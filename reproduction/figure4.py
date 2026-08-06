"""Clean-room AlphaGenome reconstruction for Figure 4B/4C: 315 bp rs12740374
core module portability into HPA liver-tissue recipient TSSs.

Ported from the working archive's report/panel_scramble_no_expression/
(run_panel_scramble_no_expression.py, run_hpa_bottom100_380bp_distance_sweep.py,
run_hpa_500_380bp_distance_sweep.py, run_hpa_liver_native_quarter.py,
run_hpa_315bp_portability_replication.py) -- none of which are part of this
repository.

Real AlphaGenome scoring throughout: every native and transfer prediction is
freshly computed (checkpointed per sequence hash so an interrupted run can
resume, but never seeded from previously-scored data).

Recipient design provenance:
- 4B ("bottom100"): the 100 recipient genes are drawn from a fixed candidate
  list hardcoded in the working archive (LOW_EXPRESSION_CANDIDATES-derived,
  not re-derivable from the archive's own source -- its original construction
  predates any traceable script) plus 3 active-liver-tissue controls
  (ALB, APOB, TTR). Their GENCODE gene/transcript/TSS resolution is
  deterministic and already recorded in
  reproduction/data/figure4b_bottom100_recipients.csv (frozen input, not an
  AlphaGenome prediction -- committed here with full provenance, same
  category as the JASPAR PFM file or Figure 3F's frozen candidate grid).
- 4C ("bottom500"/"middle500"/"top500"): fully re-derived at run time from
  the public HPA v24.1 download + GENCODE, by nTPM rank -- no frozen input.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .common import Audit, api_key, download, sha256_file
from .figure1 import _decompress, _save_svg
from .figure2 import CHROM, LIVER, RS_POS, _ag

DATA_DIR = Path(__file__).resolve().parent / "data"
BOTTOM100_RECIPIENTS_PATH = DATA_DIR / "figure4b_bottom100_recipients.csv"

HPA_URL = "https://v24.proteinatlas.org/download/rna_tissue_consensus.tsv.zip"
GENCODE_URL = "https://storage.googleapis.com/alphagenome/reference/gencode/hg38/gencode.v46.annotation.gtf.gz.feather"
# The working archive's Figure 4 scripts use this UCSC-sourced build (already
# documented as data/SOURCES.tsv's "ucsc_hg38"), not the NCBI no-alt-analysis-set
# build reproduce.py's fetch_hg38() downloads for Figures 1/3 -- confirmed to
# matter: the NCBI build has rare IUPAC ambiguity codes (Y/W/R) at a handful
# of positions that AlphaGenome's API rejects and that the UCSC build resolves
# differently, so using the wrong build is a real correctness issue here, not
# just a technical one.
UCSC_HG38_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
UCSC_HG38_EXPECTED_SHA256 = "5be01555d98347fdb3714dc84c6f77c9d8bc774adcf32c6f7a8fa06f5baf5e51"

SEQ_LEN = 2**20
TSS_HALF_WIDTH = 2_000
CORE_UPSTREAM_BP, CORE_DOWNSTREAM_BP = 179, 135
CORE_LENGTH = CORE_UPSTREAM_BP + CORE_DOWNSTREAM_BP + 1  # 315
MOTIF_START_HG38 = RS_POS - 1
MINOR_MOTIF = "GTTGCTCAAT"
DISTANCE_PANEL_DISTANCES = [0, 20, 30, 100, 1_000, 10_000, 100_000, 500_000]
COHORT_PANEL_DISTANCE = 30
COHORTS = ["bottom500", "middle500", "top500"]
DONOR_GROUPS = ["scrambled_control", "rs127_major", "rs127_minor"]


# --- shared: donor construction, sequence assembly, native/transfer states -

class GeneRecord:
    __slots__ = ("gene_symbol", "gene_id", "chrom", "strand", "start", "end", "tss",
                 "gene_type", "recipient_class", "transcript_id", "transcript_type", "tss_source", "tss_rank")

    def __init__(self, gene_symbol, gene_id, chrom, strand, start, end, tss, gene_type,
                 recipient_class, transcript_id="", transcript_type="", tss_source="gene_start", tss_rank=0):
        self.gene_symbol, self.gene_id, self.chrom, self.strand = gene_symbol, gene_id, chrom, strand
        self.start, self.end, self.tss, self.gene_type = start, end, tss, gene_type
        self.recipient_class, self.transcript_id, self.transcript_type = recipient_class, transcript_id, transcript_type
        self.tss_source, self.tss_rank = tss_source, tss_rank

    def as_dict(self) -> dict[str, object]:
        return {slot: getattr(self, slot) for slot in self.__slots__}


class DonorFragment:
    __slots__ = ("donor_id", "donor_group", "length", "sequence", "allele", "donor_start_hg38", "donor_end_hg38", "motif_status")

    def __init__(self, donor_id, donor_group, length, sequence, allele, donor_start_hg38, donor_end_hg38, motif_status):
        self.donor_id, self.donor_group, self.length = donor_id, donor_group, length
        self.sequence, self.allele = sequence, allele
        self.donor_start_hg38, self.donor_end_hg38, self.motif_status = donor_start_hg38, donor_end_hg38, motif_status


def _fig4_scramble(seq: str, *, seed: int, preserve_indices: set[int] | None = None) -> str:
    """Composition-preserving shuffle of every base in `seq` except the
    (0-based) indices in `preserve_indices`."""
    preserve_indices = preserve_indices or set()
    chars = list(seq.upper())
    idxs = [i for i, base in enumerate(chars) if i not in preserve_indices and base in "ACGT"]
    values = [chars[i] for i in idxs]
    rng = np.random.default_rng(int(seed))
    rng.shuffle(values)
    for i, value in zip(idxs, values, strict=True):
        chars[i] = value
    return "".join(chars)


def fetch_ucsc_hg38(run_dir: Path, audit: Audit) -> Path:
    raw = run_dir / "raw"
    fasta_gz = raw / "ucsc_hg38.fa.gz"
    fasta = raw / "ucsc_hg38.fa"
    with audit.step("4B/4C: download UCSC-sourced hg38 (the build the working archive actually used)"):
        info = download(UCSC_HG38_URL, fasta_gz)
        if not any(entry.get("path") == info["path"] for entry in audit.downloads):
            audit.downloads.append(info)
            audit.save()
        _decompress(fasta_gz, fasta, UCSC_HG38_EXPECTED_SHA256)
    import pysam

    if not Path(str(fasta) + ".fai").exists():
        pysam.faidx(str(fasta))
    return fasta


def _fetch_seq(fasta, chrom: str, start0: int, end0: int) -> str:
    seq = fasta.fetch(chrom, int(start0), int(end0)).upper()
    if len(seq) != int(end0) - int(start0):
        raise ValueError(f"FASTA length mismatch for {chrom}:{start0}-{end0}")
    return seq


def make_asymmetric_315_donors(fasta) -> list[DonorFragment]:
    """The three constructs displayed in Figure 4B/4C: minor (T) allele,
    major (G) allele, and a composition-matched scrambled minor control,
    each the native -179..+135 bp core around rs12740374."""
    start1 = RS_POS - CORE_UPSTREAM_BP
    end1 = RS_POS + CORE_DOWNSTREAM_BP
    major = _fetch_seq(fasta, CHROM, start1 - 1, end1)
    if len(major) != CORE_LENGTH:
        raise ValueError(f"Expected {CORE_LENGTH} bp donor, got {len(major)} bp")
    rs_idx = RS_POS - start1
    if major[rs_idx] != "G":
        raise ValueError(f"Expected G at donor index {rs_idx}, observed {major[rs_idx]}")
    minor = major[:rs_idx] + "T" + major[rs_idx + 1:]
    motif_idx = MOTIF_START_HG38 - start1
    if minor[motif_idx: motif_idx + len(MINOR_MOTIF)] != MINOR_MOTIF:
        raise ValueError("Minor donor does not contain the expected C/EBP motif")
    scrambled = _fig4_scramble(minor, seed=12740374 + CORE_LENGTH, preserve_indices={rs_idx})
    if sorted(scrambled) != sorted(minor) or scrambled[rs_idx] != "T":
        raise ValueError("Scrambled donor failed composition/allele checks")
    return [
        DonorFragment("rs127_minor_scramble_315_asym_m179_p135", "scrambled_control", CORE_LENGTH, scrambled, "minor", start1, end1, "scrambled"),
        DonorFragment("rs127_major_315_asym_m179_p135", "rs127_major", CORE_LENGTH, major, "major", start1, end1, "major_native"),
        DonorFragment("rs127_minor_315_asym_m179_p135", "rs127_minor", CORE_LENGTH, minor, "minor", start1, end1, "minor_native"),
    ]


def _safe_id(text: str) -> str:
    import re
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(text)).strip("_")


def _state_prefix(rec: GeneRecord) -> str:
    gene_id_base = rec.gene_id.split(".")[0]
    tx = rec.transcript_id.split(".")[0] if rec.transcript_id else f"rank{rec.tss_rank}"
    return _safe_id(f"{rec.gene_symbol}__{gene_id_base}__{tx}__tss{rec.tss}")


def _recipient_interval(rec: GeneRecord, seq_len: int) -> tuple[int, int]:
    start0 = int(rec.tss) - 1 - (int(seq_len) // 2)
    return start0, start0 + int(seq_len)


def _replace_segment(seq: str, *, interval_start0: int, start1: int, end1: int, replacement: str) -> str:
    idx0, idx1 = int(start1) - 1 - int(interval_start0), int(end1) - int(interval_start0)
    if idx0 < 0 or idx1 > len(seq) or idx1 <= idx0 or (idx1 - idx0) != len(replacement):
        raise ValueError(f"Replacement {start1}-{end1} does not fit interval / length mismatch.")
    return seq[:idx0] + replacement.upper() + seq[idx1:]


def _replacement_bounds_for_distance(rec: GeneRecord, length: int, distance: int) -> tuple[int, int]:
    """Upstream donor placement leaving `distance` native bp between the TSS
    and the donor's nearest edge (transcription-oriented)."""
    distance, length = int(distance), int(length)
    if rec.strand == "-":
        start1 = int(rec.tss) + distance + 1
        end1 = start1 + length - 1
    else:
        end1 = int(rec.tss) - distance - 1
        start1 = end1 - length + 1
    return start1, end1


class SequenceState:
    __slots__ = ("state_id", "gene", "sequence", "interval_start0", "interval_end0", "state_kind",
                 "donor_id", "donor_group", "donor_length", "placement", "donor_allele", "motif_status")

    def __init__(self, state_id, gene, sequence, interval_start0, interval_end0, state_kind,
                 donor_id="native", donor_group="native", donor_length=0, placement="native",
                 donor_allele="", motif_status=""):
        self.state_id, self.gene, self.sequence = state_id, gene, sequence
        self.interval_start0, self.interval_end0, self.state_kind = interval_start0, interval_end0, state_kind
        self.donor_id, self.donor_group, self.donor_length = donor_id, donor_group, donor_length
        self.placement, self.donor_allele, self.motif_status = placement, donor_allele, motif_status


def _build_native_state(fasta, rec: GeneRecord, seq_len: int) -> SequenceState:
    start0, end0 = _recipient_interval(rec, seq_len)
    seq = _fetch_seq(fasta, rec.chrom, start0, end0)
    return SequenceState(f"{_state_prefix(rec)}__native", rec, seq, start0, end0, "native")


def _build_distance_states(native: SequenceState, donors: list[DonorFragment], distances: list[int]) -> list[SequenceState]:
    rec = native.gene
    states = [native]
    for distance in distances:
        for donor in donors:
            start1, end1 = _replacement_bounds_for_distance(rec, donor.length, distance)
            seq = _replace_segment(native.sequence, interval_start0=native.interval_start0, start1=start1, end1=end1, replacement=donor.sequence)
            placement = f"upstream_{int(distance):06d}bp"
            states.append(SequenceState(
                _safe_id(f"{_state_prefix(rec)}__{placement}__{donor.donor_id}"), rec, seq,
                native.interval_start0, native.interval_end0, "transfer",
                donor_id=donor.donor_id, donor_group=donor.donor_group, donor_length=donor.length,
                placement=placement, donor_allele=donor.allele, motif_status=donor.motif_status,
            ))
    return states


def _state_metadata(state: SequenceState) -> dict[str, object]:
    rec = state.gene
    return {
        "state_id": state.state_id, "state_kind": state.state_kind, "gene_symbol": rec.gene_symbol,
        "gene_id": rec.gene_id, "gene_type": rec.gene_type, "recipient_class": rec.recipient_class,
        "transcript_id": rec.transcript_id, "transcript_type": rec.transcript_type, "tss_source": rec.tss_source,
        "tss_rank": rec.tss_rank, "chrom": rec.chrom, "strand": rec.strand, "gene_start": rec.start,
        "gene_end": rec.end, "tss_hg38": rec.tss, "interval_start0": state.interval_start0,
        "interval_end0": state.interval_end0, "donor_id": state.donor_id, "donor_group": state.donor_group,
        "donor_length": state.donor_length, "donor_allele": state.donor_allele, "motif_status": state.motif_status,
        "placement": state.placement, "sequence_sha256": hashlib.sha256(state.sequence.encode("ascii")).hexdigest(),
    }


def _mean_track_summary(track_data: Any, *, interval: Any, rec: GeneRecord, tss_half_width: int, ontology_terms: list[str]) -> dict[str, float]:
    if track_data is None:
        return {"rna_liver_all_mean": np.nan, "rna_liver_gene_strand_mean": np.nan, "rna_liver_primary": np.nan, "n_liver_tracks": 0}
    start0 = max(int(interval.start), int(rec.tss) - 1 - int(tss_half_width))
    end0 = min(int(interval.end), int(rec.tss) - 1 + int(tss_half_width) + 1)
    region = interval.__class__(chromosome=rec.chrom, start=start0, end=end0)
    sliced = track_data.slice_by_interval(region, match_resolution=True)
    if sliced is None:
        return {"rna_liver_all_mean": np.nan, "rna_liver_gene_strand_mean": np.nan, "rna_liver_primary": np.nan, "n_liver_tracks": 0}
    values = np.asarray(sliced.values, dtype=float)
    if values.ndim == 1:
        values = values[:, np.newaxis]
    meta = sliced.metadata.reset_index(drop=True)
    ont_mask = meta["ontology_curie"].astype(str).isin(ontology_terms).to_numpy() if "ontology_curie" in meta.columns else np.ones(values.shape[1], dtype=bool)
    if not np.any(ont_mask):
        ont_mask = np.ones(values.shape[1], dtype=bool)
    per_target = np.nanmean(values, axis=0)
    strand_series = meta.get("strand", pd.Series(["."] * len(meta))).fillna(".").astype(str).to_numpy()
    strand_mask = ont_mask & ((strand_series == rec.strand) | (strand_series == "."))
    all_mean = float(np.nanmean(per_target[ont_mask])) if np.any(ont_mask) else np.nan
    strand_mean = float(np.nanmean(per_target[strand_mask])) if np.any(strand_mask) else np.nan
    primary = strand_mean if np.isfinite(strand_mean) else all_mean
    return {"rna_liver_all_mean": all_mean, "rna_liver_gene_strand_mean": strand_mean, "rna_liver_primary": primary, "n_liver_tracks": int(np.sum(ont_mask))}


def _score_states(run_dir: Path, audit: Audit, panel: str, states: list[SequenceState], *, batch_size: int, max_workers: int) -> pd.DataFrame:
    """Real, fresh AlphaGenome scoring, checkpointed per (state, sequence)
    hash -- never seeded from previously-scored data."""
    genome, dna_client, dna_model, _ = _ag()
    cache = run_dir / "predictions" / f"Figure{panel}_state_cache"
    cache.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    missing: list[SequenceState] = []
    for state in states:
        digest = hashlib.sha256(state.sequence.encode("ascii")).hexdigest()
        path = cache / f"{_safe_id(state.state_id)}_{digest[:16]}.json"
        if path.exists():
            rows.append({**_state_metadata(state), **json.loads(path.read_text())})
        else:
            missing.append(state)
    client = None
    for start in range(0, len(missing), batch_size):
        batch = missing[start: start + batch_size]
        if client is None:
            client = dna_client.create(api_key(), model_version=dna_model.ModelVersion.ALL_FOLDS, timeout=300)
        with audit.step(f"{panel}: score states {start + 1}-{start + len(batch)} of {len(missing)}"):
            intervals = [genome.Interval(chromosome=s.gene.chrom, start=int(s.interval_start0), end=int(s.interval_end0)) for s in batch]
            outputs = client.predict_sequences(
                sequences=[s.sequence for s in batch], requested_outputs={dna_client.OutputType.RNA_SEQ},
                ontology_terms=[LIVER], intervals=intervals, progress_bar=False, max_workers=max_workers,
            )
            for state, interval, output in zip(batch, intervals, outputs, strict=True):
                summary = _mean_track_summary(output.rna_seq, interval=interval, rec=state.gene, tss_half_width=TSS_HALF_WIDTH, ontology_terms=[LIVER])
                digest = hashlib.sha256(state.sequence.encode("ascii")).hexdigest()
                (cache / f"{_safe_id(state.state_id)}_{digest[:16]}.json").write_text(json.dumps(summary))
                rows.append({**_state_metadata(state), **summary})
            audit.add_api_calls(panel, len(batch))
            audit.add_api_requests(panel, 1)
    return pd.DataFrame(rows)


def _add_distance_columns(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out["upstream_distance_bp"] = pd.to_numeric(out.placement.astype(str).str.extract(r"upstream_(\d+)bp", expand=False), errors="coerce")
    out.loc[out.state_kind.eq("native"), "upstream_distance_bp"] = np.nan
    return out


def _add_deltas(frame: pd.DataFrame) -> pd.DataFrame:
    pred = _add_distance_columns(frame)
    native = pred[pred.state_kind.eq("native")][["gene_symbol", "rna_liver_primary", "rna_liver_all_mean", "rna_liver_gene_strand_mean"]].rename(
        columns={"rna_liver_primary": "native_rna_liver_primary", "rna_liver_all_mean": "native_rna_liver_all_mean", "rna_liver_gene_strand_mean": "native_rna_liver_gene_strand_mean"}
    )
    pred = pred.merge(native, on="gene_symbol", how="left")
    pred["delta_liver_primary_vs_native"] = pred.rna_liver_primary - pred.native_rna_liver_primary
    return pred


# --- 4B: bottom-100 recipient design (frozen) + eight-distance sweep -------

def _read_bottom100_recipients() -> list[GeneRecord]:
    design = pd.read_csv(BOTTOM100_RECIPIENTS_PATH)
    return [
        GeneRecord(str(r.gene_symbol), str(r.gene_id), str(r.chrom), str(r.strand), int(r.start), int(r.end),
                   int(r.tss), str(r.gene_type), str(r.recipient_class), str(r.transcript_id), str(r.transcript_type),
                   str(r.tss_source), int(r.tss_rank))
        for r in design.itertuples(index=False)
    ]


def run_fig4b(run_dir: Path, audit: Audit, *, batch_size: int = 128, max_workers: int = 8) -> None:
    import pysam

    ucsc_fasta_path = fetch_ucsc_hg38(run_dir, audit)
    with audit.step("4B: build 315bp donor + bottom100 distance-sweep design"):
        fasta = pysam.FastaFile(str(ucsc_fasta_path))
        donors = make_asymmetric_315_donors(fasta)
        recipients = [r for r in _read_bottom100_recipients() if r.recipient_class == "hpa_liver_bottom100"]
        states: list[SequenceState] = []
        for rec in recipients:
            native = _build_native_state(fasta, rec, SEQ_LEN)
            states.extend(_build_distance_states(native, donors, DISTANCE_PANEL_DISTANCES))
        out = run_dir / "derived/Figure4B_distance_response"
        out.mkdir(parents=True, exist_ok=True)
        pd.DataFrame([_state_metadata(s) for s in states]).to_csv(out / "sequence_state_design.csv", index=False)

    scored = _score_states(run_dir, audit, "4B", states, batch_size=batch_size, max_workers=max_workers)
    with audit.step("4B: compute deltas and summarize by distance"):
        pred = _add_deltas(scored)
        pred.to_csv(out / "predictions_with_deltas.csv", index=False)
        transfer = pred[pred.state_kind.eq("transfer") & pred.donor_group.isin(["rs127_major", "rs127_minor"])]
        summary_rows = []
        for (distance, group_name), group in transfer.groupby(["upstream_distance_bp", "donor_group"]):
            vals = group.delta_liver_primary_vs_native.replace([np.inf, -np.inf], np.nan).dropna()
            summary_rows.append({
                "upstream_distance_bp": float(distance), "donor_group": group_name, "n": int(len(vals)),
                "mean": float(vals.mean()) if len(vals) else np.nan, "median": float(vals.median()) if len(vals) else np.nan,
                "sem": float(vals.sem()) if len(vals) > 1 else np.nan,
            })
        summary = pd.DataFrame(summary_rows).sort_values(["upstream_distance_bp", "donor_group"]).reset_index(drop=True)
        summary.to_csv(out / "Figure4B_distance_response.tsv", sep="\t", index=False)
        _render4b(summary, run_dir / "figures/Figure4B.svg")


def _render4b(summary: pd.DataFrame, output: Path) -> None:
    colors = {"rs127_major": "#8A99A8", "rs127_minor": "#3F7FE5"}
    fig, ax = plt.subplots(figsize=(3.4, 2.4))
    for group, color in colors.items():
        sub = summary[summary.donor_group.eq(group)].sort_values("upstream_distance_bp")
        ax.errorbar(sub.upstream_distance_bp.clip(lower=0.5), sub["mean"], yerr=sub["sem"], marker="o", markersize=3, color=color, label=group)
    ax.set_xscale("symlog", linthresh=1)
    ax.set_xlabel("Upstream distance from TSS (bp)")
    ax.set_ylabel("Δ liver RNA vs native")
    ax.legend(frameon=False, fontsize=7)
    _save_svg(fig, output)


# --- 4C: HPA bottom/middle/top-500 cohorts (re-derived from public HPA) ---

HPA_EXPECTED_SHA256 = "cdedaeaf3cdfc89e22b3891ea24ae2afabc0afd26d8883076121a363608450b6"


def _fetch_hpa(run_dir: Path, audit: Audit, hpa_file: Path | None = None) -> Path:
    """The Human Protein Atlas download endpoint blocks scripted requests
    (HTTP 403 regardless of User-Agent/Referer -- confirmed 2026-08-05, same
    class of issue as the already-documented Wang 2018 supplement). Same
    fallback convention as --wang-xls: accept a manually downloaded copy,
    checksum-verified against the value recorded in data/SOURCES.tsv."""
    raw = run_dir / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    tsv_path = raw / "rna_tissue_consensus.tsv"
    if tsv_path.exists():
        return tsv_path
    if hpa_file is not None:
        import shutil

        from .common import sha256_file

        observed = sha256_file(hpa_file)
        if observed != HPA_EXPECTED_SHA256:
            raise ValueError(f"--hpa-file checksum mismatch: expected {HPA_EXPECTED_SHA256}, got {observed}")
        with audit.step("4C: stage manually supplied HPA v24.1 consensus tissue RNA expression"):
            shutil.copy2(hpa_file, tsv_path)
            info = {"url": f"supplied:{hpa_file}", "path": str(tsv_path), "bytes": tsv_path.stat().st_size, "sha256": observed, "reused": False}
            audit.downloads.append(info)
            audit.save()
        return tsv_path
    zip_path = raw / "rna_tissue_consensus.tsv.zip"
    with audit.step("4C: download HPA v24.1 consensus tissue RNA expression"):
        info = download(HPA_URL, zip_path)
        if not any(entry.get("path") == info["path"] for entry in audit.downloads):
            audit.downloads.append(info)
            audit.save()
    import zipfile
    with zipfile.ZipFile(zip_path) as zf:
        names = zf.namelist()
        with zf.open(names[0]) as src, tsv_path.open("wb") as dst:
            dst.write(src.read())
    return tsv_path


_STANDARD_CHROMOSOMES = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def _load_hpa_liver(path: Path, gene_rows: pd.DataFrame) -> pd.DataFrame:
    hpa = pd.read_csv(path, sep="\t")
    liver = hpa[hpa["Tissue"].astype(str).str.lower().eq("liver")].copy()
    liver["hpa_gene_id_base"] = liver["Gene"].astype(str).str.split(".").str[0]
    liver["hpa_gene_name_sort"] = liver["Gene name"].astype(str).str.upper()
    liver["hpa_liver_ntpm"] = pd.to_numeric(liver["nTPM"], errors="coerce")
    liver = liver.dropna(subset=["hpa_liver_ntpm"])
    liver = liver.sort_values(["hpa_gene_id_base", "hpa_gene_name_sort"]).drop_duplicates("hpa_gene_id_base", keep="first").reset_index(drop=True)
    # Restrict to genes GENCODE resolves on a standard chromosome *before*
    # ranking -- otherwise unresolvable/non-standard-chromosome genes shift
    # the bottom/middle/top-500 cutoffs by a handful of genes at the edges.
    eligible = set(gene_rows[gene_rows.Chromosome.astype(str).isin(_STANDARD_CHROMOSOMES)].gene_id_base.astype(str))
    liver = liver[liver.hpa_gene_id_base.astype(str).isin(eligible)].copy()
    return liver.sort_values(["hpa_liver_ntpm", "hpa_gene_name_sort", "hpa_gene_id_base"]).reset_index(drop=True)


def _select_hpa_cohort(hpa_liver: pd.DataFrame, cohort: str, size: int = 500) -> pd.DataFrame:
    if cohort == "bottom500":
        selected = hpa_liver.head(size).copy()
    elif cohort == "top500":
        selected = hpa_liver.tail(size).iloc[::-1].copy()
    elif cohort == "middle500":
        start = max(0, (len(hpa_liver) - size) // 2)
        selected = hpa_liver.iloc[start: start + size].copy()
    else:
        raise ValueError(f"Unsupported cohort: {cohort}")
    return selected.reset_index(drop=True)


def _ensure_gencode(run_dir: Path, audit: Audit) -> Path:
    path = run_dir / "raw" / "gencode.v46.annotation.gtf.gz.feather"
    if path.exists():
        return path
    with audit.step("4C: download GENCODE v46 annotation feather"):
        info = download(GENCODE_URL, path)
        if not any(entry.get("path") == info["path"] for entry in audit.downloads):
            audit.downloads.append(info)
            audit.save()
    return path


def _load_gencode_tables(gencode_feather: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    cols = ["Chromosome", "Feature", "Start", "End", "Strand", "gene_name", "gene_id", "gene_type", "transcript_id", "transcript_type"]
    gencode = pd.read_feather(gencode_feather, columns=cols)
    gencode["gene_id_base"] = gencode.gene_id.astype(str).str.split(".").str[0]
    gene_rows = gencode[gencode.Feature.eq("gene")].copy()
    gene_rows["is_pc"] = gene_rows.gene_type.astype(str).eq("protein_coding")
    gene_rows = gene_rows.sort_values(["is_pc", "Chromosome", "Start"], ascending=[False, True, True]).drop_duplicates("gene_id_base", keep="first")
    tx_rows = gencode[gencode.Feature.eq("transcript")].copy()
    return gene_rows, tx_rows


def _cohort_recipients(hpa_selected: pd.DataFrame, gene_rows: pd.DataFrame, tx_rows: pd.DataFrame, *, recipient_class: str) -> list[GeneRecord]:
    merged = hpa_selected.merge(gene_rows, left_on="hpa_gene_id_base", right_on="gene_id_base", how="left")
    merged = merged[merged.gene_id.notna()].copy()
    tx = tx_rows[tx_rows.gene_id_base.isin(set(merged.hpa_gene_id_base.astype(str)))].copy()
    tx["candidate_tss"] = np.where(tx.Strand.astype(str).eq("+"), tx.Start, tx.End).astype(int)
    tx["is_pc_tx"] = tx.transcript_type.astype(str).eq("protein_coding")
    tx = tx.sort_values(["gene_id_base", "is_pc_tx", "Start", "End"], ascending=[True, False, True, True])
    records: list[GeneRecord] = []
    for _, row in merged.iterrows():
        gene_id_base = str(row.hpa_gene_id_base)
        candidates = tx[tx.gene_id_base.eq(gene_id_base)]
        protein_coding = candidates[candidates.is_pc_tx]
        chosen = protein_coding if not protein_coding.empty else candidates
        if chosen.empty:
            continue
        first = chosen.iloc[0]
        records.append(GeneRecord(
            str(row.gene_name).upper(), str(row.gene_id), str(row.Chromosome), str(row.Strand),
            int(row.Start), int(row.End), int(first.candidate_tss), str(row.gene_type), recipient_class,
            str(first.transcript_id), str(first.transcript_type),
            "protein_coding_transcript_tss" if str(first.transcript_type) == "protein_coding" else "transcript_tss", 1,
        ))
    return records


def run_fig4c(run_dir: Path, audit: Audit, *, batch_size: int = 128, max_workers: int = 8, hpa_file: Path | None = None) -> None:
    import pysam

    ucsc_fasta_path = fetch_ucsc_hg38(run_dir, audit)
    hpa_path = _fetch_hpa(run_dir, audit, hpa_file)
    gencode_path = _ensure_gencode(run_dir, audit)
    with audit.step("4C: select HPA bottom/middle/top-500 cohorts and build 30bp donor-transfer design"):
        gene_rows, tx_rows = _load_gencode_tables(gencode_path)
        hpa_liver = _load_hpa_liver(hpa_path, gene_rows)
        fasta = pysam.FastaFile(str(ucsc_fasta_path))
        donors = make_asymmetric_315_donors(fasta)
        out = run_dir / "derived/Figure4C_foldchange_cohorts"
        out.mkdir(parents=True, exist_ok=True)
        cohort_states: dict[str, list[SequenceState]] = {}
        skipped: list[dict[str, object]] = []
        for cohort in COHORTS:
            selected = _select_hpa_cohort(hpa_liver, cohort)
            recipients = _cohort_recipients(selected, gene_rows, tx_rows, recipient_class=f"hpa_liver_{cohort}")
            states: list[SequenceState] = []
            for rec in recipients:
                # A TSS near a chromosome end can put the +/-524288bp scoring
                # window outside the chromosome (e.g. chr12 in hg38 is
                # 133,275,309 bp; a subtelomeric TSS's window can overrun
                # that by several hundred kb). The working archive skips such
                # genes rather than failing the run; matched here.
                try:
                    native = _build_native_state(fasta, rec, SEQ_LEN)
                except ValueError as exc:
                    skipped.append({"cohort": cohort, "gene_symbol": rec.gene_symbol, "chrom": rec.chrom, "tss": rec.tss, "reason": str(exc)})
                    continue
                states.extend(_build_distance_states(native, donors, [COHORT_PANEL_DISTANCE]))
            cohort_states[cohort] = states
        if skipped:
            pd.DataFrame(skipped).to_csv(out / "skipped_recipients.csv", index=False)
        all_states = [s for states in cohort_states.values() for s in states]

    scored = _score_states(run_dir, audit, "4C", all_states, batch_size=batch_size, max_workers=max_workers)
    with audit.step("4C: compute fold changes by cohort"):
        cohort_by_gene = {}
        for cohort, states in cohort_states.items():
            for s in states:
                cohort_by_gene[s.gene.gene_symbol] = cohort
        scored["cohort"] = scored.gene_symbol.map(cohort_by_gene)
        pred = _add_deltas(scored)
        transfer = pred[pred.state_kind.eq("transfer") & pred.donor_group.isin(DONOR_GROUPS)].copy()
        transfer = transfer[(transfer.rna_liver_primary > 0) & (transfer.native_rna_liver_primary > 0)].copy()
        transfer["fold_change_vs_native"] = transfer.rna_liver_primary / transfer.native_rna_liver_primary
        transfer["log2_fold_change_vs_native"] = np.log2(transfer.fold_change_vs_native)
        transfer.to_csv(out / "Figure4C_fold_change_values.csv", index=False)
        summary_rows = []
        for (cohort, donor_group), group in transfer.groupby(["cohort", "donor_group"]):
            values = group.log2_fold_change_vs_native
            summary_rows.append({
                "cohort": cohort, "donor_group": donor_group, "n": int(len(values)),
                "mean": float(values.mean()), "median": float(values.median()),
                "q25": float(values.quantile(0.25)), "q75": float(values.quantile(0.75)),
                "minimum": float(values.min()), "maximum": float(values.max()),
                "fraction_up": float((values > 0).mean()), "fraction_down": float((values < 0).mean()),
            })
        summary = pd.DataFrame(summary_rows).sort_values(["cohort", "donor_group"]).reset_index(drop=True)
        summary.to_csv(out / "Figure4C_fold_change_summary.csv", index=False)
        _render4c(summary, run_dir / "figures/Figure4C.svg")


def _render4c(summary: pd.DataFrame, output: Path) -> None:
    colors = {"rs127_minor": "#3F7FE5", "rs127_major": "#8A99A8", "scrambled_control": "#F77F2F"}
    fig, ax = plt.subplots(figsize=(2.8, 2.05))
    x = np.arange(len(COHORTS))
    width = 0.25
    for i, group in enumerate(DONOR_GROUPS):
        sub = summary[summary.donor_group.eq(group)].set_index("cohort").reindex(COHORTS)
        ax.bar(x + (i - 1) * width, sub["mean"], width, color=colors[group], label=group)
    ax.set_xticks(x, COHORTS, rotation=0, fontsize=7)
    ax.set_ylabel("log2 fold change vs native")
    ax.axhline(0, color="#333333", linewidth=0.7)
    ax.legend(frameon=False, fontsize=6)
    _save_svg(fig, output)
