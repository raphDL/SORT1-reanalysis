"""Clean-room AlphaGenome pipelines for Figure 2B, 2C, 2E, and 2F."""

from __future__ import annotations

import hashlib
import json
import math
import shutil
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from scipy import stats

from .common import Audit, api_key, download, make_alphagenome_importable, sha256_file, utc_now
from .figure1 import EXPECTED_UNCOMPRESSED, HG38_URL, _decompress


CHROM = "chr1"
RS_POS = 109_274_968
RS_MAJOR = "G"
RS_MINOR = "T"
SEQ_LEN = 524_288
LIVER = "UBERON:0002107"
HEPG2 = "EFO:0001187"
CUT_LEFT = 109_274_966
MOTIF_START = 109_274_967
MOTIF_END = 109_274_976
TSS_HALF_WIDTH = 2_000
GENES = ("SORT1", "PSRC1", "CELSR2")
BASES = ("A", "C", "G", "T")

WANG_URL = (
    "https://www.ahajournals.org/action/downloadSupplement?"
    "doi=10.1161%2FATVBAHA.117.310103&file=atv310103_ds.xls"
)
WANG_SHA256 = "99a2d4218e3e3a9afdd04ecded3a38878c001272f42eaadc99b8c3caeeefa049"
KIRCHER_COMMIT = "05d2ffb965090d3f5dd27dfb038cec493a15ab35"
KIRCHER_URL = (
    "https://raw.githubusercontent.com/kircherlab/MPRA_SaturationMutagenesis/"
    f"{KIRCHER_COMMIT}/data/elements.tsv.gz"
)
KIRCHER_SHA256 = "fec2eed91fe27af3aae07ebce2eca65e9bad4bb6abba5d8c27f478887dd7b134"

GENE_TABLE = pd.DataFrame(
    [
        ("SORT1", "SORT1", "ENSG00000134243.12", 109_397_918, "-"),
        ("PSRC1", "PSRC1", "ENSG00000134222.16", 109_283_186, "-"),
        ("CELSR2", "CELSR2", "ENSG00000143126.8", 109_249_538, "+"),
        ("SARS", "SARS1", "ENSG00000031698.13", 109_213_917, "+"),
        ("PSMA5", "PSMA5", "ENSG00000143106.13", 109_426_448, "-"),
    ],
    columns=["gene_symbol", "gencode_gene_name", "gene_id", "tss_hg38", "gene_strand"],
)


def _ag() -> tuple[Any, Any, Any, Any]:
    make_alphagenome_importable()
    from alphagenome.data import genome
    from alphagenome.models import dna_client, dna_model, variant_scorers

    return genome, dna_client, dna_model, variant_scorers


def fetch_hg38(run_dir: Path, audit: Audit) -> Path:
    raw = run_dir / "raw"
    fasta_gz = raw / "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
    item = download(HG38_URL, fasta_gz)
    if not any(entry.get("path") == item["path"] for entry in audit.downloads):
        audit.downloads.append(item)
        audit.save()
    fasta = raw / "hg38.fa"
    _decompress(fasta_gz, fasta, EXPECTED_UNCOMPRESSED["hg38.fa"])
    import pysam

    if not Path(str(fasta) + ".fai").exists():
        pysam.faidx(str(fasta))
    return fasta


def fetch_kircher(run_dir: Path, audit: Audit) -> Path:
    destination = run_dir / "raw" / "kircher_elements.tsv.gz"
    item = download(KIRCHER_URL, destination)
    if item["sha256"] != KIRCHER_SHA256:
        raise ValueError("The downloaded Kircher elements file does not match the pinned checksum")
    if not any(entry.get("path") == item["path"] for entry in audit.downloads):
        audit.downloads.append(item)
        audit.save()
    return destination


def stage_wang(run_dir: Path, audit: Audit, supplied: Path | None = None) -> Path:
    """Stage the original Wang spreadsheet, allowing a manual publisher download."""
    destination = run_dir / "raw" / "atv310103_ds.xls"
    if destination.exists():
        observed = sha256_file(destination)
        reused = True
        source = (
            f"manual_original_file:{supplied.expanduser().resolve().name}"
            if supplied is not None else WANG_URL
        )
    elif supplied is not None:
        supplied = supplied.expanduser().resolve()
        if not supplied.exists():
            raise FileNotFoundError(f"Wang spreadsheet not found: {supplied}")
        shutil.copy2(supplied, destination)
        observed = sha256_file(destination)
        reused = False
        source = f"manual_original_file:{supplied.name}"
    else:
        try:
            item = download(WANG_URL, destination)
            observed = str(item["sha256"])
            reused = bool(item["reused"])
            source = WANG_URL
        except RuntimeError as exc:
            raise RuntimeError(
                "The Wang publisher download is bot-gated. Download the original "
                "atv310103_ds.xls from DOI 10.1161/ATVBAHA.117.310103 and rerun with "
                "--wang-xls /path/to/atv310103_ds.xls; its SHA-256 will be verified."
            ) from exc
    if observed != WANG_SHA256:
        destination.unlink(missing_ok=True)
        raise ValueError(f"Wang spreadsheet checksum mismatch: expected {WANG_SHA256}, got {observed}")
    item = {
        "url": source,
        "publisher_url": WANG_URL,
        "path": str(destination),
        "bytes": destination.stat().st_size,
        "sha256": observed,
        "reused": reused,
    }
    if not any(entry.get("path") == item["path"] for entry in audit.downloads):
        audit.downloads.append(item)
        audit.save()
    return destination


def prepare_figure2_inputs(
    run_dir: Path, audit: Audit, panels: list[str], wang_xls: Path | None = None
) -> dict[str, Path]:
    result: dict[str, Path] = {}
    if {"2B", "2C"} & set(panels):
        result["fasta"] = fetch_hg38(run_dir, audit)
    if "2B" in panels:
        result["wang"] = stage_wang(run_dir, audit, wang_xls)
    if {"2E", "2F"} & set(panels):
        result["kircher"] = fetch_kircher(run_dir, audit)
    return result


def _mean_rna(track: Any, interval: Any, genes: pd.DataFrame = GENE_TABLE) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for gene in genes.itertuples(index=False):
        region = interval.__class__(
            chromosome=CHROM,
            start=max(int(interval.start), int(gene.tss_hg38) - TSS_HALF_WIDTH),
            end=min(int(interval.end), int(gene.tss_hg38) + TSS_HALF_WIDTH + 1),
        )
        sliced = track.slice_by_interval(region, match_resolution=True)
        values = np.asarray(sliced.values, dtype=float)
        if values.ndim == 1:
            values = values[:, None]
        metadata = sliced.metadata.copy()
        metadata["target_index"] = np.arange(len(metadata))
        metadata["gene_symbol"] = gene.gene_symbol
        metadata["gencode_gene_name"] = gene.gencode_gene_name
        metadata["gene_id"] = gene.gene_id
        metadata["gene_tss_hg38"] = int(gene.tss_hg38)
        metadata["tss_half_width"] = TSS_HALF_WIDTH
        metadata["rna_mean_tss_pm2kb"] = np.nanmean(values, axis=0)
        rows.append(metadata)
    return pd.concat(rows, ignore_index=True)


@dataclass(frozen=True)
class SequenceState:
    state_id: str
    sequence: str
    sequence_sha256: str


def _score_sequences(
    states: list[SequenceState], run_dir: Path, audit: Audit, panel: str,
    *, batch_size: int, max_workers: int,
) -> pd.DataFrame:
    genome, dna_client, dna_model, _ = _ag()
    start0 = (RS_POS - 1) - SEQ_LEN // 2
    interval = genome.Interval(CHROM, start0, start0 + SEQ_LEN)
    cache = run_dir / "predictions" / f"Figure{panel}_sequence_cache"
    cache.mkdir(parents=True, exist_ok=True)
    frames: list[pd.DataFrame] = []
    missing: list[SequenceState] = []
    for state in states:
        path = cache / f"{state.sequence_sha256}.tsv"
        if path.exists():
            frames.append(pd.read_csv(path, sep="\t"))
        else:
            missing.append(state)
    client = None
    for offset in range(0, len(missing), batch_size):
        batch = missing[offset : offset + batch_size]
        if client is None:
            client = dna_client.create(
                api_key(), model_version=dna_model.ModelVersion.ALL_FOLDS, timeout=300
            )
        outputs = client.predict_sequences(
            sequences=[state.sequence for state in batch],
            requested_outputs={dna_client.OutputType.RNA_SEQ},
            ontology_terms=[LIVER],
            intervals=[interval] * len(batch),
            progress_bar=False,
            max_workers=max_workers,
        )
        audit.add_api_requests(panel, 1)
        audit.add_api_calls(panel, len(batch))
        for state, output in zip(batch, outputs, strict=True):
            frame = _mean_rna(output.rna_seq, interval)
            frame["state_id"] = state.state_id
            frame["sequence_sha256"] = state.sequence_sha256
            frame.to_csv(cache / f"{state.sequence_sha256}.tsv", sep="\t", index=False)
            frames.append(frame)
    if not frames:
        raise RuntimeError(f"No AlphaGenome RNA predictions were generated for Figure {panel}")
    return pd.concat(frames, ignore_index=True)


def _replace_rs(sequence: str, index: int, allele: str) -> str:
    if sequence[index].upper() != RS_MAJOR:
        raise ValueError(f"Expected {RS_MAJOR} at rs12740374, found {sequence[index]}")
    return sequence[:index] + allele + sequence[index + 1 :]


def _parse_wang(path: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    excel = pd.ExcelFile(path)
    for sheet in excel.sheet_names:
        table = pd.read_excel(path, sheet_name=sheet)
        table.columns = [str(value).strip() for value in table.columns]
        if "Position (hg18)" not in table or "Cons" not in table:
            raise ValueError(f"Unexpected columns in Wang sheet {sheet}")
        for column in ("Reads", "VarFreq", "IndelFreq", "Position (hg18)"):
            table[column] = pd.to_numeric(table[column], errors="coerce")
        modeled = table[table.Cons.astype(str).str.match(r"^\*/[+-][ACGT]+$", na=False)].copy()
        modeled = modeled.sort_values("Reads", ascending=False).reset_index(drop=True)
        denominator = float(modeled.Reads.sum())
        experiment = sheet.replace("CRISPR-SNP ", "").replace(" ", "_")
        for rank, row in enumerate(modeled.itertuples(index=False), start=1):
            cons = str(row.Cons)
            operation = "deletion" if cons[2] == "-" else "insertion"
            event = cons[3:].upper()
            # Namedtuple field normalization varies across pandas versions; use the table row below.
            source = modeled.iloc[rank - 1]
            position = int(source["Position (hg18)"])
            edit_id = f"{experiment}_rank{rank:03d}_{position}_{cons.replace('*/', '')}"
            rows.append(
                {
                    "experiment": experiment,
                    "rank_in_sheet": rank,
                    "edit_id": edit_id.replace("+", "_").replace("/", "_"),
                    "position_hg18": position,
                    "reference_base_reported": str(source.Ref).upper(),
                    "reads": int(source.Reads),
                    "VarFreq": float(source.VarFreq),
                    "IndelFreq": float(source.IndelFreq),
                    "Cons": cons,
                    "operation": operation,
                    "event_sequence": event,
                    "event_length": len(event),
                    "p_i": float(source.Reads) / denominator,
                }
            )
    return pd.DataFrame(rows)


def _wang_sequences(indels: pd.DataFrame, fasta_path: Path) -> tuple[pd.DataFrame, dict[str, str]]:
    import pysam

    fasta = pysam.FastaFile(str(fasta_path))
    start0 = (RS_POS - 1) - SEQ_LEN // 2
    extra = int(indels.event_length.max()) + 64
    reference = fasta.fetch(CHROM, start0, start0 + SEQ_LEN + extra).upper()
    rs_index = RS_POS - 1 - start0
    minor = _replace_rs(reference, rs_index, RS_MINOR)
    offset = RS_POS - 109_619_113
    hashes: list[str] = []
    sequences: dict[str, str] = {}
    matches: list[bool] = []
    for row in indels.itertuples(index=False):
        anchor = int(row.position_hg18) + offset
        anchor_index = anchor - 1 - start0
        boundary = anchor_index + 1
        if row.operation == "deletion":
            observed = minor[boundary : boundary + int(row.event_length)]
            edited = minor[:boundary] + minor[boundary + int(row.event_length) :]
            matches.append(observed == row.event_sequence)
        else:
            edited = minor[:boundary] + row.event_sequence + minor[boundary:]
            matches.append(True)
        sequence = edited[:SEQ_LEN]
        digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
        hashes.append(digest)
        sequences[digest] = sequence
    out = indels.copy()
    out["sequence_sha256"] = hashes
    out["event_sequence_matches_minor_background"] = matches
    if not out.loc[out.operation.eq("deletion"), "event_sequence_matches_minor_background"].all():
        raise ValueError("At least one Wang deletion does not match the downloaded GRCh38 sequence")
    major_seq = reference[:SEQ_LEN]
    minor_seq = minor[:SEQ_LEN]
    sequences[hashlib.sha256(major_seq.encode("ascii")).hexdigest()] = major_seq
    sequences[hashlib.sha256(minor_seq.encode("ascii")).hexdigest()] = minor_seq
    return out, sequences


def _labels(rows: pd.DataFrame) -> list[str]:
    seen: dict[str, int] = {}
    labels: list[str] = []
    for row in rows.itertuples(index=False):
        label = str(row.Cons)
        if len(label) > 14:
            label = f"-{int(row.event_length)}bp" if row.operation == "deletion" else f"+{int(row.event_length)}bp"
        seen[label] = seen.get(label, 0) + 1
        labels.append(label if seen[label] == 1 else f"{label} ({seen[label]})")
    return labels


def run_fig2b(
    run_dir: Path, audit: Audit, fasta: Path, wang: Path,
    *, batch_size: int = 8, max_workers: int = 4,
) -> None:
    indels, sequences = _wang_sequences(_parse_wang(wang), fasta)
    chosen = pd.concat(
        [
            indels[indels.experiment.eq(experiment)].sort_values("rank_in_sheet").head(50)
            for experiment in ("human_hepatocytes", "mouse")
        ],
        ignore_index=True,
    )
    chosen["experiment_label"] = chosen.experiment.map(
        {"human_hepatocytes": "Human", "mouse": "Mouse"}
    )
    columns = (
        chosen.groupby("sequence_sha256", sort=False)
        .agg(
            best_rank=("rank_in_sheet", "min"),
            max_p_i=("p_i", "max"),
            total_reads=("reads", "sum"),
            Cons=("Cons", "first"),
            operation=("operation", "first"),
            event_sequence=("event_sequence", "first"),
            event_length=("event_length", "first"),
            edit_ids=("edit_id", lambda values: ";".join(values.astype(str))),
            source_experiments=("experiment_label", lambda values: ",".join(sorted(set(values.astype(str))))),
        )
        .reset_index()
        .sort_values(["best_rank", "max_p_i", "Cons"], ascending=[True, False, True], kind="mergesort")
        .head(50)
        .reset_index(drop=True)
    )
    if len(columns) != 50:
        raise ValueError(f"Expected 50 unique Wang repair sequences, found {len(columns)}")
    columns["fused_rank"] = np.arange(1, len(columns) + 1)
    columns["column_label"] = _labels(columns)
    # Locate baselines by sequence construction rather than any publication table.
    import pysam

    reference = pysam.FastaFile(str(fasta)).fetch(
        CHROM, (RS_POS - 1) - SEQ_LEN // 2, (RS_POS - 1) - SEQ_LEN // 2 + SEQ_LEN
    ).upper()
    rs_index = SEQ_LEN // 2
    minor_sequence = _replace_rs(reference, rs_index, RS_MINOR)
    major_hash = hashlib.sha256(reference.encode("ascii")).hexdigest()
    minor_hash = hashlib.sha256(minor_sequence.encode("ascii")).hexdigest()
    states = [SequenceState("major", sequences[major_hash], major_hash), SequenceState("minor", sequences[minor_hash], minor_hash)]
    states.extend(
        SequenceState(f"edit_{digest[:16]}", sequences[digest], digest)
        for digest in columns.sequence_sha256.astype(str)
    )
    rna = _score_sequences(states, run_dir, audit, "2B", batch_size=batch_size, max_workers=max_workers)
    target = rna[rna.target_index.astype(int).eq(4) & rna.gene_symbol.isin(GENES)].copy()
    baseline = target[target.sequence_sha256.eq(minor_hash)].set_index("gene_symbol").rna_mean_tss_pm2kb
    target["percent_change"] = [
        100.0 * (value - baseline[gene]) / baseline[gene]
        for gene, value in zip(target.gene_symbol, target.rna_mean_tss_pm2kb, strict=False)
    ]
    matrix = (
        target[target.sequence_sha256.isin(columns.sequence_sha256)]
        .pivot_table(index="gene_symbol", columns="sequence_sha256", values="percent_change", aggfunc="mean")
        .reindex(index=GENES, columns=columns.sequence_sha256)
    )
    matrix.columns = columns.column_label
    if matrix.isna().any().any():
        raise ValueError("Incomplete Figure 2B matrix")
    output = run_dir / "derived" / "Figure2B_top50_repair_outcomes"
    output.mkdir(parents=True, exist_ok=True)
    matrix.reset_index().to_csv(output / "matrix.csv", index=False)
    columns.to_csv(output / "columns.csv", index=False)
    indels.to_csv(output / "wang_indels_reconstructed.tsv", sep="\t", index=False)
    _render_heatmap(matrix, run_dir / "figures" / "Figure2B.svg", "Top repair outcomes")


def _build_deletion_grid(fasta_path: Path) -> tuple[pd.DataFrame, list[SequenceState], str]:
    import pysam

    fasta = pysam.FastaFile(str(fasta_path))
    start0 = (RS_POS - 1) - SEQ_LEN // 2
    extended = fasta.fetch(CHROM, start0, start0 + SEQ_LEN + 84).upper()
    rs_index = RS_POS - 1 - start0
    cut_index = CUT_LEFT - start0
    minor = _replace_rs(extended, rs_index, RS_MINOR)
    major = extended[:SEQ_LEN]
    minor_seq = minor[:SEQ_LEN]
    rows: list[dict[str, object]] = []
    sequences: dict[str, str] = {}

    def add(design_id: str, operation: str, sequence: str, up: int, down: int) -> None:
        digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
        sequences[digest] = sequence
        rows.append(
            {
                "design_id": design_id,
                "operation": operation,
                "length": up + down,
                "upstream_bases": float(up),
                "downstream_bases": float(down),
                "geometry": "full_xy_grid" if operation == "deletion" else "baseline",
                "sequence_sha256": digest,
            }
        )

    add("major", "baseline_major", major, 0, 0)
    add("minor", "baseline_minor", minor_seq, 0, 0)
    for up in range(11):
        for down in range(11):
            edited = minor[: cut_index - up] + minor[cut_index + down :]
            add(f"del_xy_u{up:02d}_d{down:02d}", "deletion", edited[:SEQ_LEN], up, down)
    designs = pd.DataFrame(rows)
    states = [
        SequenceState(f"seq_{digest[:16]}", sequence, digest)
        for digest, sequence in sorted(sequences.items())
    ]
    return designs, states, hashlib.sha256(minor_seq.encode("ascii")).hexdigest()


def _minor_baseline(expanded: pd.DataFrame, keys: list[str], minor_hash: str) -> pd.DataFrame:
    """Build a one-row-per-`keys` baseline table for the minor-allele sequence.

    `expanded` can legitimately contain more than one design row for the same
    underlying sequence -- e.g. the (upstream=0, downstream=0) grid cell
    deletes zero bases and so hashes identically to the standalone "minor"
    design. `.drop_duplicates(keys)` alone assumes every row sharing a `keys`
    combination is byte-identical, which does not hold if an optional
    metadata column is NaN for some rows (pandas groups NaN keys together,
    so a `validate="many_to_one"` merge downstream can intermittently raise
    "Merge keys are not unique in right dataset"). Aggregate explicitly
    instead, which is unique by construction, and fail loudly with a
    diagnosable message if the underlying values genuinely disagree rather
    than raising an opaque pandas MergeError.
    """
    subset = expanded[expanded.sequence_sha256.eq(minor_hash)][keys + ["rna_mean_tss_pm2kb"]]
    grouped = subset.groupby(keys, as_index=False, dropna=False)
    spread = grouped["rna_mean_tss_pm2kb"].agg(lambda values: float(values.max() - values.min()))
    inconsistent = spread[spread.rna_mean_tss_pm2kb > 1e-9]
    if not inconsistent.empty:
        raise ValueError(
            "Minor-allele baseline disagrees across duplicate design rows for the same "
            f"{keys} combination (max spread {inconsistent.rna_mean_tss_pm2kb.max():.3g}); "
            f"first offending row: {inconsistent.iloc[0].to_dict()}"
        )
    return grouped["rna_mean_tss_pm2kb"].mean().rename(
        columns={"rna_mean_tss_pm2kb": "rna_mean_tss_pm2kb_minor"}
    )


def run_fig2c(
    run_dir: Path, audit: Audit, fasta: Path, *, batch_size: int = 8, max_workers: int = 4
) -> None:
    designs, states, minor_hash = _build_deletion_grid(fasta)
    rna = _score_sequences(states, run_dir, audit, "2C", batch_size=batch_size, max_workers=max_workers)
    expanded = designs.merge(rna, on="sequence_sha256", how="left", validate="many_to_many")
    keys = [
        column for column in (
            "gene_symbol", "gencode_gene_name", "target_index", "name", "strand",
            "ontology_curie", "biosample_name", "Assay title", "tss_half_width",
        ) if column in expanded
    ]
    baseline = _minor_baseline(expanded, keys, minor_hash)
    expanded = expanded.merge(baseline, on=keys, how="left", validate="many_to_one")
    expanded["delta_vs_minor"] = expanded.rna_mean_tss_pm2kb - expanded.rna_mean_tss_pm2kb_minor
    expanded["percent_change_vs_minor"] = expanded.delta_vs_minor / expanded.rna_mean_tss_pm2kb_minor.replace(0, np.nan)
    output = run_dir / "derived" / "Figure2C_deletion_grid.csv"
    expanded.to_csv(output, index=False)
    _render_deletion_grid(expanded, run_dir / "figures" / "Figure2C.svg")


def _centered(values: np.ndarray) -> TwoSlopeNorm:
    finite = values[np.isfinite(values)]
    maximum = float(np.max(np.abs(finite))) if finite.size else 1.0
    maximum = maximum or 1.0
    return TwoSlopeNorm(vmin=-maximum, vcenter=0, vmax=maximum)


def _render_heatmap(matrix: pd.DataFrame, path: Path, title: str) -> None:
    values = matrix.to_numpy(float)
    fig, ax = plt.subplots(figsize=(max(8.7, 0.18 * matrix.shape[1] + 1.9), 2.65))
    image = ax.pcolormesh(values, cmap="RdBu_r", norm=_centered(values), edgecolors="white", linewidth=0.35)
    ax.set_ylim(matrix.shape[0], 0)
    ax.set_title(title)
    ax.set_yticks(np.arange(matrix.shape[0]) + 0.5, matrix.index)
    ax.set_xticks(np.arange(matrix.shape[1]) + 0.5, matrix.columns, rotation=90, fontsize=6)
    fig.colorbar(image, ax=ax, label="RNA change vs minor (%)")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def _render_deletion_grid(data: pd.DataFrame, path: Path) -> None:
    subset = data[
        data.geometry.eq("full_xy_grid") & data.target_index.astype(int).eq(4) & data.gene_symbol.isin(GENES)
    ].copy()
    subset["display"] = 100.0 * subset.percent_change_vs_minor
    fig, axes = plt.subplots(1, 3, figsize=(10.6, 3.7), constrained_layout=True)
    image = None
    for ax, gene in zip(axes, GENES, strict=True):
        matrix = subset[subset.gene_symbol.eq(gene)].pivot_table(
            index="upstream_bases", columns="downstream_bases", values="display", aggfunc="mean"
        ).reindex(index=range(11), columns=range(11))
        image = ax.imshow(matrix, origin="lower", cmap="RdBu_r", norm=_centered(subset.display.to_numpy()))
        ax.set_title(gene)
        ax.set_xlabel("bases deleted downstream")
    axes[0].set_ylabel("bases deleted upstream")
    fig.colorbar(image, ax=axes, label="RNA change vs minor (%)")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def _construct_alleles(data: pd.DataFrame) -> pd.DataFrame:
    out = data.copy()
    out["GenomeRef"] = out.Ref.astype(str).str.upper()
    out["GenomeAlt"] = out.Alt.astype(str).str.upper()
    out["ConstructRef"] = [RS_MINOR if int(pos) == RS_POS else ref for pos, ref in zip(out.Position, out.GenomeRef)]
    out["ConstructAlt"] = [
        RS_MAJOR if int(pos) == RS_POS and alt == RS_MINOR else alt
        for pos, alt in zip(out.Position, out.GenomeAlt)
    ]
    out["variant_id_genome"] = out.Position.astype(int).astype(str) + ":" + out.GenomeRef + ">" + out.GenomeAlt
    out["variant_id_construct"] = out.Position.astype(int).astype(str) + ":" + out.ConstructRef + ">" + out.ConstructAlt
    out["is_rs12740374_position"] = out.Position.astype(int).eq(RS_POS)
    out["is_rs12740374_construct_to_genome"] = out.is_rs12740374_position & out.ConstructAlt.eq(RS_MAJOR)
    return out


def load_sort1_mpra(path: Path) -> pd.DataFrame:
    raw = pd.read_csv(path, sep="\t", low_memory=False).rename(
        columns={"Chrom": "Chromosome", "Pos": "Position", "Barcodes": "Tags", "Coefficient": "Value", "pValue": "P-Value"}
    )
    raw = raw[
        raw.Release.eq("GRCh38")
        & raw.Element.isin(("SORT1", "SORT1.2", "SORT1-flip"))
        & raw.Ref.astype(str).str.upper().isin(BASES)
        & raw.Alt.astype(str).str.upper().isin(BASES)
    ].copy()
    raw.Position = raw.Position.astype(int)
    # The study's released SORT1-specific table reports the MPRA coefficient to
    # two decimals; reconstruct that public representation from elements.tsv.gz.
    raw["Value"] = pd.to_numeric(raw["Value"], errors="coerce").round(2)
    raw = _construct_alleles(raw)
    raw["offset"] = raw.Position - RS_POS
    primary = raw[raw.Element.isin(("SORT1", "SORT1.2"))]
    consensus = (
        primary.groupby(
            ["Position", "GenomeRef", "GenomeAlt", "ConstructRef", "ConstructAlt", "variant_id_genome",
             "variant_id_construct", "offset", "is_rs12740374_position", "is_rs12740374_construct_to_genome"],
            as_index=False,
        )
        .agg(
            kircher_primary_log2_effect=("Value", "mean"), kircher_primary_sd=("Value", "std"),
            kircher_primary_n_constructs=("Value", "size"), kircher_primary_min_p=("P-Value", "min"),
            kircher_primary_tags=("Tags", "sum"), kircher_primary_dna=("DNA", "sum"), kircher_primary_rna=("RNA", "sum"),
        )
        .fillna({"kircher_primary_sd": 0.0})
        .sort_values(["Position", "ConstructAlt"])
        .reset_index(drop=True)
    )
    consensus["kircher_primary_abs_log2_effect"] = consensus.kircher_primary_log2_effect.abs()
    if len(consensus) != 1798:
        raise ValueError(f"Expected 1,798 SORT1 MPRA substitutions, found {len(consensus)}")
    return consensus


def _position_chunks(positions: list[int], size: int) -> list[tuple[int, int]]:
    chunks: list[tuple[int, int]] = []
    run: list[int] = []
    previous = None
    for position in positions:
        if previous is None or (position == previous + 1 and len(run) < size):
            run.append(position)
        else:
            chunks.append((run[0], run[-1] + 1))
            run = [position]
        previous = position
    if run:
        chunks.append((run[0], run[-1] + 1))
    return chunks


def _variant_record(variant: Any) -> dict[str, object] | None:
    if variant is None:
        return None
    position = int(variant.position)
    ref = str(variant.reference_bases).upper()
    alt = str(variant.alternate_bases).upper()
    if len(ref) != 1 or len(alt) != 1 or ref == alt or ref not in BASES or alt not in BASES:
        return None
    genome_ref = RS_MAJOR if position == RS_POS else ref
    genome_alt = RS_MINOR if position == RS_POS and alt == RS_MAJOR else alt
    return {
        "Position": position, "GenomeRef": genome_ref, "GenomeAlt": genome_alt,
        "ConstructRef": ref, "ConstructAlt": alt,
        "variant_id_genome": f"{position}:{genome_ref}>{genome_alt}",
        "variant_id_construct": f"{position}:{ref}>{alt}",
        "is_rs12740374_position": position == RS_POS,
        "is_rs12740374_construct_to_genome": position == RS_POS and alt == RS_MAJOR,
        "ag_api_ref": ref, "ag_api_alt": alt,
    }


def _score_sort1_modality(
    consensus: pd.DataFrame, run_dir: Path, audit: Audit, modality: str,
    *, chunk_size: int, max_workers: int,
    model_version_name: str = "ALL_FOLDS", panel: str = "2E", cache_dirname: str = "Figure2E_kircher",
) -> pd.DataFrame:
    genome, dna_client, dna_model, variant_scorers = _ag()
    cache = run_dir / "predictions" / cache_dirname / "cache"
    cache.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    positions = sorted(consensus.Position.unique())
    center = (positions[0] + positions[-1]) // 2
    window = 131_072 if modality == "rna" else 16_384
    interval = genome.Interval(CHROM, center, center).resize(window)
    scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS[
        {"rna": "RNA_SEQ", "atac": "ATAC", "h3k27ac": "CHIP_HISTONE"}[modality]
    ]
    interval_variant = genome.Variant(CHROM, RS_POS, RS_MAJOR, RS_MINOR, "rs12740374_T_background")
    client = dna_client.create(api_key(), model_version=getattr(dna_model.ModelVersion, model_version_name), timeout=300)
    for start, end in _position_chunks(positions, chunk_size):
        path = cache / f"{modality}_{start}_{end - 1}.tsv"
        paths.append(path)
        if path.exists():
            continue
        for attempt in range(1, 6):
            try:
                outputs = client.score_ism_variants(
                    interval=interval,
                    ism_interval=genome.Interval(CHROM, start - 1, end - 1),
                    interval_variant=interval_variant,
                    variant_scorers=[scorer],
                    progress_bar=False,
                    max_workers=max_workers,
                )
                break
            except Exception:
                if attempt == 5:
                    raise
                time.sleep(min(2**attempt, 20))
        audit.add_api_requests(panel, 1)
        rows: list[dict[str, object]] = []
        for variant_outputs in outputs:
            if not variant_outputs:
                continue
            adata = variant_outputs[0]
            record = _variant_record(adata.uns.get("variant"))
            if record is None:
                continue
            if modality == "rna":
                genes = adata.obs[adata.obs.gene_name.isin(GENES)]
                tracks = adata.var[adata.var.ontology_curie.astype(str).eq(HEPG2)]
                for gene_index, gene in genes.iterrows():
                    strand = str(gene.get("strand", ""))
                    for track_index, track in tracks.iterrows():
                        track_strand = str(track.get("strand", ""))
                        if strand in {"+", "-"} and track_strand not in {strand, "."}:
                            continue
                        rows.append(
                            {**record, "gene_name": str(gene.gene_name), "gene_id": str(gene.get("gene_id", "")),
                             "ontology_curie": HEPG2, "biosample_name": str(track.get("biosample_name", "")),
                             "assay": str(track.get("Assay title", track.get("assay", ""))),
                             "track_strand": track_strand, "raw_score": float(adata[gene_index, track_index].X[0, 0])}
                        )
            else:
                metadata = adata.var.copy()
                mask = metadata.ontology_curie.astype(str).eq(HEPG2).to_numpy()
                if modality == "h3k27ac":
                    mark = next(
                        (column for column in ("histone_mark_code", "histone_mark") if column in metadata),
                        None,
                    )
                    mask &= (
                        metadata[mark].astype(str).eq("H3K27ac").to_numpy()
                        if mark is not None else np.zeros(len(metadata), dtype=bool)
                    )
                values = np.asarray(adata.X, dtype=float).reshape(1, -1)[0, mask]
                if values.size == 0:
                    continue
                prefix = f"ag_{modality}"
                rows.append(
                    {**record, f"{prefix}_n_tracks": int(values.size),
                     f"{prefix}_mean_score": float(np.nanmean(values)),
                     f"{prefix}_abs_score": float(np.nanmean(np.abs(values))),
                     f"{prefix}_max_abs_score": float(np.nanmax(np.abs(values)))}
                )
        pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
        audit.add_api_calls(panel, len({row["variant_id_construct"] for row in rows}))
    return pd.concat([pd.read_csv(path, sep="\t") for path in paths], ignore_index=True)


def _aggregate_sort1_rna(raw: pd.DataFrame) -> pd.DataFrame:
    keys = ["Position", "GenomeRef", "GenomeAlt", "ConstructRef", "ConstructAlt", "variant_id_genome", "variant_id_construct"]
    per_gene = raw.groupby(keys + ["gene_name"], as_index=False).agg(
        ag_rna_lnfc=("raw_score", "mean"), ag_rna_n_tracks=("raw_score", "size")
    )
    per_gene["ag_rna_percent_change"] = 100.0 * np.expm1(per_gene.ag_rna_lnfc)
    out = per_gene.pivot_table(index=keys, columns="gene_name", values=["ag_rna_lnfc", "ag_rna_n_tracks", "ag_rna_percent_change"]).reset_index()
    out.columns = ["_".join(str(part) for part in col if str(part)).strip("_") if isinstance(col, tuple) else str(col) for col in out.columns]
    for gene in GENES:
        out[f"ag_rna_abs_percent_change_{gene}"] = out[f"ag_rna_percent_change_{gene}"].abs()
    signed = [f"ag_rna_percent_change_{gene}" for gene in GENES]
    absolute = [f"ag_rna_abs_percent_change_{gene}" for gene in GENES]
    out["ag_rna_3gene_mean_percent_change"] = out[signed].mean(axis=1)
    out["ag_rna_3gene_mean_abs_percent_change"] = out[absolute].mean(axis=1)
    out["ag_rna_3gene_max_abs_percent_change"] = out[absolute].max(axis=1)
    return out


def run_fig2e(
    run_dir: Path, audit: Audit, kircher: Path, *, chunk_size: int = 10, max_workers: int = 1
) -> None:
    consensus = load_sort1_mpra(kircher)
    frames = {
        modality: _score_sort1_modality(consensus, run_dir, audit, modality, chunk_size=chunk_size, max_workers=max_workers)
        for modality in ("atac", "h3k27ac", "rna")
    }
    keys = ["Position", "GenomeRef", "GenomeAlt", "ConstructRef", "ConstructAlt", "variant_id_genome", "variant_id_construct"]
    matched = consensus.copy()
    for frame in (frames["atac"], frames["h3k27ac"], _aggregate_sort1_rna(frames["rna"])):
        extra = [column for column in frame if column not in keys and column not in matched]
        matched = matched.merge(frame[keys + extra].drop_duplicates(keys), on=keys, how="left")
    if len(matched.dropna(subset=["ag_atac_mean_score", "ag_h3k27ac_mean_score", "ag_rna_3gene_mean_percent_change"])) != 1798:
        raise ValueError("Figure 2E did not produce 1,798 complete matched substitutions")
    output = run_dir / "derived" / "Figure2E_kircher_correlations.tsv"
    matched.to_csv(output, sep="\t", index=False)
    _render_fig2e(matched, run_dir / "figures" / "Figure2E.svg")


def _render_fig2e(data: pd.DataFrame, path: Path) -> None:
    panels = [
        ("RNA", "ag_rna_3gene_mean_percent_change"),
        ("ATAC", "ag_atac_mean_score"),
        ("H3K27ac", "ag_h3k27ac_mean_score"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(6.5, 2.2))
    for ax, (title, column) in zip(axes, panels, strict=True):
        motif = data.Position.between(MOTIF_START, MOTIF_END)
        ax.scatter(data.loc[~motif, "kircher_primary_log2_effect"], data.loc[~motif, column], s=3, alpha=.25)
        ax.scatter(data.loc[motif, "kircher_primary_log2_effect"], data.loc[motif, column], s=5, color="#c93636")
        x = data.kircher_primary_log2_effect.to_numpy(float)
        y = data[column].to_numpy(float)
        slope, intercept = np.polyfit(x, y, 1)
        grid = np.linspace(x.min(), x.max(), 100)
        ax.plot(grid, slope * grid + intercept, color="black", lw=.8)
        ax.set_title(f"{title}  r={stats.pearsonr(x, y).statistic:.2f}")
        ax.set_xlabel("Kircher log2 effect")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


ELEMENTS = {
    "F9": (("F9",), HEPG2, ("DNASE", "ATAC"), "HepG2"),
    "FOXE1": (("FOXE1",), "EFO:0002791", ("DNASE",), "HeLa-S3"),
    "LDLR": (("LDLR", "LDLR.2"), HEPG2, ("DNASE", "ATAC"), "HepG2"),
    "MYC_rs11986220": (("MYCrs11986220",), "EFO:0005726", ("DNASE",), "LNCaP clone FGC"),
    "PKLR": (("PKLR-48h",), "EFO:0002067", ("DNASE", "ATAC"), "K562"),
    "SORT1": (("SORT1", "SORT1.2"), HEPG2, ("DNASE", "ATAC"), "HepG2"),
}


def load_multielement(path: Path) -> dict[str, pd.DataFrame]:
    raw = pd.read_csv(path, sep="\t", low_memory=False)
    raw = raw[
        raw.Release.eq("GRCh38") & raw.Ref.astype(str).str.len().eq(1)
        & raw.Alt.astype(str).str.len().eq(1) & (pd.to_numeric(raw.Barcodes, errors="coerce") >= 10)
    ].copy()
    raw.Pos = raw.Pos.astype(int)
    raw.Chrom = "chr" + raw.Chrom.astype(str).str.replace(r"\.0$", "", regex=True)
    result = {}
    for element, (experiments, _, _, _) in ELEMENTS.items():
        subset = raw[raw.Element.isin(experiments)]
        consensus = subset.groupby(["Chrom", "Pos", "Ref", "Alt"], as_index=False).agg(
            mpra_effect=("Coefficient", "mean"), mpra_sd=("Coefficient", "std"),
            mpra_n_experiments=("Coefficient", "size"), mpra_min_barcodes=("Barcodes", "min"),
            mpra_min_p=("pValue", "min"),
        ).fillna({"mpra_sd": 0.0})
        consensus["element"] = element
        result[element] = consensus
    return result


def _score_multielement(
    client: Any, element: str, experimental: pd.DataFrame, modality: str,
    run_dir: Path, audit: Audit,
) -> pd.DataFrame:
    genome, dna_client, _, variant_scorers = _ag()
    ontology = ELEMENTS[element][1]
    cache = run_dir / "predictions" / "Figure2F_kircher_multielement" / modality.lower()
    cache.mkdir(parents=True, exist_ok=True)
    chrom = str(experimental.Chrom.iloc[0])
    start, end = int(experimental.Pos.min()), int(experimental.Pos.max())
    interval = genome.Interval(chrom, (start + end) // 2, (start + end) // 2).resize(16_384)
    output_type = dna_client.OutputType.ATAC if modality == "ATAC" else dna_client.OutputType.DNASE
    scorer = variant_scorers.CenterMaskScorer(output_type, None, variant_scorers.AggregationType.DIFF_LOG2_SUM)
    paths = []
    for chunk_start in range(start, end + 1, 10):
        chunk_end = min(chunk_start + 10, end + 1)
        path = cache / f"{element}_{chunk_start}_{chunk_end - 1}.tsv"
        paths.append(path)
        if path.exists():
            continue
        for attempt in range(1, 6):
            try:
                outputs = client.score_ism_variants(
                    interval=interval, ism_interval=genome.Interval(chrom, chunk_start, chunk_end),
                    variant_scorers=[scorer], progress_bar=False, max_workers=1,
                )
                break
            except Exception:
                if attempt == 5:
                    raise
                time.sleep(min(2**attempt, 20))
        audit.add_api_requests("2F", 1)
        rows = []
        for variant_outputs in outputs:
            adata = variant_outputs[0]
            variant = adata.uns["variant"]
            mask = adata.var.ontology_curie.astype(str).eq(ontology).to_numpy()
            values = np.asarray(adata.X, dtype=float).reshape(1, -1)[0, mask]
            rows.append(
                {"Chrom": chrom, "Pos": int(variant.position), "Ref": str(variant.reference_bases).upper(),
                 "Alt": str(variant.alternate_bases).upper(),
                 "ag_accessibility_effect": float(np.nanmean(values)), "n_matched_tracks": int(values.size),
                 "ontology": ontology, "modality": modality}
            )
        pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
        audit.add_api_calls("2F", len(rows))
    return pd.concat([pd.read_csv(path, sep="\t") for path in paths], ignore_index=True).drop_duplicates(
        ["Chrom", "Pos", "Ref", "Alt"]
    )


def _random_effects(table: pd.DataFrame) -> dict[str, float]:
    z = np.arctanh(table.spearman_rho.clip(-.999999, .999999).to_numpy())
    variance = 1.0 / (table.n.to_numpy(float) - 3.0)
    fixed = 1.0 / variance
    mean = float(np.sum(fixed * z) / np.sum(fixed))
    q = float(np.sum(fixed * (z - mean) ** 2))
    df = len(z) - 1
    c = float(np.sum(fixed) - np.sum(fixed**2) / np.sum(fixed))
    tau2 = max(0.0, (q - df) / c)
    weights = 1.0 / (variance + tau2)
    pooled = float(np.sum(weights * z) / np.sum(weights))
    se = float(np.sqrt(1.0 / np.sum(weights)))
    return {
        "k_elements": len(z), "pooled_spearman_rho": float(np.tanh(pooled)),
        "pooled_95ci_low": float(np.tanh(pooled - 1.96 * se)),
        "pooled_95ci_high": float(np.tanh(pooled + 1.96 * se)), "tau2_fisher_z": tau2,
        "Q": q, "Q_df": df, "Q_p": float(stats.chi2.sf(q, df)),
        "I2_percent": float(max(0.0, (q - df) / q) * 100.0) if q else 0.0,
    }


def _analyze_multielement(merged: pd.DataFrame, kircher: Path) -> tuple[pd.DataFrame, dict[str, Any]]:
    rows = []
    for (element, modality), group in merged.groupby(["element", "modality"], sort=False):
        clean = group[["mpra_effect", "ag_accessibility_effect"]].dropna()
        rho = stats.spearmanr(clean.mpra_effect, clean.ag_accessibility_effect)
        pearson = stats.pearsonr(clean.mpra_effect, clean.ag_accessibility_effect)
        z = np.arctanh(np.clip(float(rho.statistic), -.999999, .999999))
        se = 1 / np.sqrt(len(clean) - 3)
        rows.append(
            {"element": element, "n": len(clean), "spearman_rho": float(rho.statistic),
             "spearman_p": float(rho.pvalue), "spearman_ci_low": float(np.tanh(z - 1.96 * se)),
             "spearman_ci_high": float(np.tanh(z + 1.96 * se)), "pearson_r": float(pearson.statistic),
             "pearson_p": float(pearson.pvalue),
             "direction_agreement": float((np.sign(clean.mpra_effect) == np.sign(clean.ag_accessibility_effect)).mean()),
             "modality": modality, "cell_line": ELEMENTS[element][3]}
        )
    table = pd.DataFrame(rows)
    by_modality = {modality: _random_effects(group) for modality, group in table.groupby("modality")}
    atac_elements = set(table.loc[table.modality.eq("ATAC"), "element"])
    matched = {
        modality: _random_effects(table[table.modality.eq(modality) & table.element.isin(atac_elements)])
        for modality in ("DNASE", "ATAC")
    }
    summary = {
        "data_sha256": sha256_file(kircher), "unique_elements": int(table.element.nunique()),
        "element_modality_pairs": len(table), "substitutions_total": int(table.n.sum()),
        "random_effects_by_modality": by_modality,
        "matched_four_random_effects_by_modality": matched,
        "matched_four_elements": sorted(atac_elements),
        "median_element_spearman": float(table.spearman_rho.median()),
        "positive_element_correlations": int((table.spearman_rho > 0).sum()),
        "analysis_mode": "all-folds", "model_version_by_element": {element: "ALL_FOLDS" for element in ELEMENTS},
    }
    return table, summary


def run_fig2f(run_dir: Path, audit: Audit, kircher: Path) -> None:
    _, dna_client, dna_model, _ = _ag()
    experimental = load_multielement(kircher)
    parts = []
    for element, table in experimental.items():
        client = dna_client.create(api_key(), model_version=dna_model.ModelVersion.ALL_FOLDS, timeout=300)
        for modality in ELEMENTS[element][2]:
            scored = _score_multielement(client, element, table, modality, run_dir, audit)
            parts.append(table.merge(scored, on=["Chrom", "Pos", "Ref", "Alt"], how="inner"))
    merged = pd.concat(parts, ignore_index=True)
    statistics, summary = _analyze_multielement(merged, kircher)
    output = run_dir / "derived" / "Figure2F_kircher_multielement"
    output.mkdir(parents=True, exist_ok=True)
    merged.to_csv(output / "kircher_multielement_matched_scores.tsv", sep="\t", index=False)
    statistics.to_csv(output / "kircher_multielement_element_statistics.tsv", sep="\t", index=False)
    (output / "RESULTS.json").write_text(json.dumps(summary, indent=2) + "\n")
    _render_fig2f(statistics, summary, run_dir / "figures" / "Figure2F.svg")


def _render_fig2f(table: pd.DataFrame, summary: dict[str, Any], path: Path) -> None:
    colors = {"ATAC": "#1b9e77", "DNASE": "#7570b3"}
    markers = {"HepG2": "o", "K562": "s", "HeLa-S3": "^", "LNCaP clone FGC": "v"}
    order = table[table.modality.eq("DNASE")].sort_values("spearman_rho").element.tolist()
    ymap = {element: index + 2 for index, element in enumerate(order)}
    fig, ax = plt.subplots(figsize=(5.8, 3.25))
    for row in table.itertuples(index=False):
        y = ymap[row.element] + {"DNASE": -.11, "ATAC": .11}[row.modality]
        ax.hlines(y, row.spearman_ci_low, row.spearman_ci_high, color=colors[row.modality])
        ax.scatter(row.spearman_rho, y, color=colors[row.modality], marker=markers[row.cell_line])
    for modality, offset in (("DNASE", -.09), ("ATAC", .09)):
        pooled = summary["matched_four_random_effects_by_modality"][modality]
        ax.hlines(.65 + offset, pooled["pooled_95ci_low"], pooled["pooled_95ci_high"], color=colors[modality])
        ax.scatter(pooled["pooled_spearman_rho"], .65 + offset, color=colors[modality], marker="D")
    ax.set_yticks([.65] + [ymap[element] for element in order], ["Matched subset (k=4)"] + order)
    ax.axvline(0, color="gray", ls="--")
    ax.set_xlabel("Spearman rho: MPRA vs predicted accessibility")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
