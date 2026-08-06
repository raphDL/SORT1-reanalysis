"""Clean-room AlphaGenome reconstruction for Figure 4H: exhaustive +/-50kb
single-nucleotide mutagenesis scan around rs12740374, scored for coordinated
SORT1/PSRC1/CELSR2 RNA(TSS) upregulation in liver, CD8+ memory T cells, and
CD14+ monocytes.

Ported from the working archive's shared ISM engine
(sort1_figure_2e_100kb_rna_ism.py::run_stage2_snv_scan / _candidate_variants_
from_windows / _mean_track_signal / _compute_gene_tss_signals, invoked
exhaustively over the full region by results/figure2/panel_a/
run_panel_a_full_region_exhaustive_snv.py for liver, and by
results/figure3/panelB/run_panel_b_figure2A_true_liver_tcell_exhaustive.py /
results/figure3/panelD's macrophage counterpart for the other two tissues),
plus the "best alt per position" synergy reduction from
results/figure2/panel_a/run_panel_a_full_region_snv_synergy.py. None of
these are part of this repository.

Exhaustive: every position in the +/-50kb window x each of its 3 non-
reference bases (300,003 variants) is scored independently per tissue via
predict_variants -- no coarse-to-fine search like Figure 3A. This is a
materially larger computation than every other Figure 4 panel done so far
(~900,000 AlphaGenome variant scores total); every one of them is freshly
computed, checkpointed by TSV-append per completed batch so an interrupted
multi-day run can resume without rescoring anything already done, but never
seeded from previously-scored data.

The reference table (17.3MB) is Zenodo-pending -- not committed to this
repository's git history (see outputs/run_manifests/
zenodo_pending_large_outputs.tsv) -- so compare_fig4h in report.py reads it
from the working archive's own copy after verifying its checksum matches
the recorded one, the same file either way.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .common import Audit, api_key
from .figure2 import CHROM, GENES, RS_POS, LIVER, _ag
from .figure3 import _gene_tss_table
from .figure4 import fetch_ucsc_hg38

REGION_HALF_WIDTH = 50_000
REGION_START = RS_POS - REGION_HALF_WIDTH
REGION_END = RS_POS + REGION_HALF_WIDTH
SEQ_WINDOW = 2 ** 19
TSS_HALF_WIDTH = 2_000
BASES = ("A", "C", "G", "T")

# Tissue tracks: liver averages every AlphaGenome track under its ontology
# (matching the archive's liver-exhaustive engine); T cell and CD14+
# monocyte pin one specific named track each (matching the archive's
# T-cell/macrophage exhaustive scripts' strict selection), since those
# ontologies were disambiguated by name in the working archive.
TISSUES: dict[str, dict[str, Any]] = {
    "liver": {"label": "liver", "ontology_curie": LIVER, "track_name": None},
    "cd14_monocyte": {"label": "CD14+ monocyte", "ontology_curie": "CL:0001054", "track_name": "CL:0001054 polyA plus RNA-seq"},
    "tcell": {"label": "T cell", "ontology_curie": "CL:0000909", "track_name": "CL:0000909 polyA plus RNA-seq"},
}
TRACK_COLORS = {"liver": "#d62728", "cd14_monocyte": "#f2b701", "tcell": "#1f78b4"}


def _candidate_variants(genome_mod, fasta, *, chrom: str, region_start: int, region_end: int) -> list[tuple[int, str, str]]:
    """(position, ref_base, alt_base) for every position in [region_start,
    region_end] (1-based, inclusive) x its 3 non-reference bases."""
    seq = fasta.fetch(chrom, region_start - 1, region_end).upper()
    if len(seq) != region_end - region_start + 1:
        raise ValueError(f"FASTA length mismatch for {chrom}:{region_start}-{region_end}")
    variants: list[tuple[int, str, str]] = []
    for offset, ref in enumerate(seq):
        if ref not in BASES:
            continue
        position = region_start + offset
        for alt in BASES:
            if alt != ref:
                variants.append((position, ref, alt))
    return variants


def _track_signal(track_data, *, genome_mod, chrom: str, center_pos: int, half_width: int, interval, ontology_curie: str, track_name: str | None) -> float:
    if track_data is None:
        return float("nan")
    region_start = max(int(interval.start), center_pos - half_width)
    region_end = min(int(interval.end), center_pos + half_width + 1)
    if region_end <= region_start:
        return float("nan")
    region = genome_mod.Interval(chromosome=chrom, start=region_start, end=region_end)
    sliced = track_data.slice_by_interval(region, match_resolution=True)
    if sliced is None:
        return float("nan")
    values = np.asarray(sliced.values, dtype=float)
    if values.ndim == 1:
        values = values[:, np.newaxis]
    if values.size == 0:
        return float("nan")
    meta = sliced.metadata
    if track_name is not None:
        mask = meta["name"].astype(str).eq(track_name).to_numpy() if "name" in meta.columns else np.zeros(values.shape[1], dtype=bool)
    else:
        mask = meta["ontology_curie"].astype(str).eq(ontology_curie).to_numpy() if "ontology_curie" in meta.columns else np.zeros(values.shape[1], dtype=bool)
    if not np.any(mask):
        return float("nan")
    return float(np.nanmean(values[:, mask]))


def _append_rows_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    frame = pd.DataFrame(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not path.exists()
    frame.to_csv(path, sep="\t", index=False, mode="a", header=write_header)


def _exhaustive_ism_scan(
    run_dir: Path, audit: Audit, *, tissue_key: str, variants: list[tuple[int, str, str]],
    gene_tss: dict[str, int], out_tsv: Path, batch_size: int, max_workers: int,
) -> pd.DataFrame:
    genome_mod, dna_client, dna_model, _ = _ag()
    spec = TISSUES[tissue_key]
    interval = genome_mod.Interval(chromosome=CHROM, start=RS_POS, end=RS_POS).resize(SEQ_WINDOW)

    done_keys: set[str] = set()
    if out_tsv.exists():
        old = pd.read_csv(out_tsv, sep="\t", usecols=["variant_key"])
        if not old.empty:
            done_keys = set(old["variant_key"].astype(str).unique().tolist())
    pending = [(pos, ref, alt) for pos, ref, alt in variants if f"{CHROM}:{pos}:{ref}>{alt}" not in done_keys]

    if pending:
        client = dna_client.create(api_key(), model_version=dna_model.ModelVersion.ALL_FOLDS, timeout=180)
        n_total = len(variants)
        for start in range(0, len(pending), batch_size):
            chunk = pending[start: start + batch_size]
            chunk_variants = [
                genome_mod.Variant(chromosome=CHROM, position=pos, reference_bases=ref, alternate_bases=alt, name=f"{CHROM}:{pos}:{ref}>{alt}")
                for pos, ref, alt in chunk
            ]
            with audit.step(f"4H {tissue_key}: score variants {start + 1}-{start + len(chunk)} of {len(pending)} (of {n_total} total)"):
                outputs = client.predict_variants(
                    intervals=interval, variants=chunk_variants, requested_outputs={dna_client.OutputType.RNA_SEQ},
                    ontology_terms=[spec["ontology_curie"]], progress_bar=False, max_workers=max_workers,
                )
                rows: list[dict[str, Any]] = []
                for (pos, ref, alt), output in zip(chunk, outputs, strict=True):
                    for gene, tss in gene_tss.items():
                        ref_val = _track_signal(
                            output.reference.rna_seq, genome_mod=genome_mod, chrom=CHROM, center_pos=int(tss),
                            half_width=TSS_HALF_WIDTH, interval=interval, ontology_curie=spec["ontology_curie"], track_name=spec["track_name"],
                        )
                        alt_val = _track_signal(
                            output.alternate.rna_seq, genome_mod=genome_mod, chrom=CHROM, center_pos=int(tss),
                            half_width=TSS_HALF_WIDTH, interval=interval, ontology_curie=spec["ontology_curie"], track_name=spec["track_name"],
                        )
                        rows.append({
                            "variant_key": f"{CHROM}:{pos}:{ref}>{alt}", "position": pos, "ref_base": ref, "alt_base": alt,
                            "gene": gene, "ref_rna_tss": ref_val, "alt_rna_tss": alt_val, "delta_rna_tss": float(alt_val - ref_val),
                        })
                _append_rows_tsv(out_tsv, rows)
                audit.add_api_calls("4H", len(chunk))
                audit.add_api_requests("4H", 1)

    return pd.read_csv(out_tsv, sep="\t")


def _best_alt_by_position(long_df: pd.DataFrame, genes: tuple[str, ...]) -> pd.DataFrame:
    wide = long_df.pivot_table(index=["variant_key", "position", "ref_base", "alt_base"], columns="gene", values="delta_rna_tss", aggfunc="mean").reset_index()
    wide.columns.name = None
    arr = wide[list(genes)].to_numpy(dtype=float)
    all_positive = np.all(arr > 0.0, axis=1)
    min_delta = np.min(arr, axis=1)
    mean_delta = np.mean(arr, axis=1)
    wide["min_delta"] = min_delta
    wide["mean_delta"] = mean_delta
    wide["synergy_score"] = np.where(all_positive, min_delta, 0.0)
    ranked = wide.sort_values(["position", "synergy_score", "min_delta", "mean_delta"], ascending=[True, False, False, False])
    return ranked.groupby("position", as_index=False).head(1).sort_values("position").reset_index(drop=True)


def run_fig4h(run_dir: Path, audit: Audit, *, batch_size: int = 32, max_workers: int = 8) -> None:
    import pysam

    ucsc_fasta_path = fetch_ucsc_hg38(run_dir, audit)
    fasta = pysam.FastaFile(str(ucsc_fasta_path))
    genome_mod, *_ = _ag()

    with audit.step("4H: build the +/-50kb x 3-alt candidate-variant universe"):
        variants = _candidate_variants(genome_mod, fasta, chrom=CHROM, region_start=REGION_START, region_end=REGION_END)
        gene_tss = {gene: tss for gene, tss in _gene_tss_table().items() if gene in GENES}

    out = run_dir / "derived/Figure4H_regional_tissue_scan"
    out.mkdir(parents=True, exist_ok=True)
    predictions_dir = run_dir / "predictions/Figure4H_ism"
    predictions_dir.mkdir(parents=True, exist_ok=True)

    combined_rows: list[pd.DataFrame] = []
    for tissue_key, spec in TISSUES.items():
        raw_tsv = predictions_dir / f"{tissue_key}_ism_tss_deltas.tsv"
        long_df = _exhaustive_ism_scan(
            run_dir, audit, tissue_key=tissue_key, variants=variants, gene_tss=gene_tss,
            out_tsv=raw_tsv, batch_size=batch_size, max_workers=max_workers,
        )
        with audit.step(f"4H {tissue_key}: best-alt-per-position synergy reduction"):
            best = _best_alt_by_position(long_df, GENES)
            best.to_csv(out / f"{tissue_key}_best_alt_by_position.tsv", sep="\t", index=False)
            track = best[["position", "synergy_score"]].copy()
            track["track"] = tissue_key
            track["track_label"] = spec["label"]
            track["position_offset_bp"] = track["position"] - RS_POS
            track["position_offset_kb"] = track["position_offset_bp"] / 1000.0
            track["synergy_x1e3"] = track["synergy_score"] * 1000.0
            combined_rows.append(track)

    with audit.step("4H: assemble combined liver/CD14-monocyte/T-cell synergy track table"):
        combined = pd.concat(combined_rows, ignore_index=True).sort_values(["track", "position"]).reset_index(drop=True)
        combined = combined[["position", "synergy_score", "track", "track_label", "position_offset_bp", "position_offset_kb", "synergy_x1e3"]]
        combined.to_csv(out / "Figure4H_regional_tissue_scan.tsv", sep="\t", index=False)
        _render4h(combined, run_dir / "figures/Figure4H.svg")


def _render4h(combined: pd.DataFrame, output: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from .figure1 import _save_svg

    fig, ax = plt.subplots(figsize=(5.0, 2.2))
    for tissue_key, spec in TISSUES.items():
        sub = combined[combined.track.eq(tissue_key)].sort_values("position")
        ax.plot(sub.position_offset_kb, sub.synergy_x1e3, color=TRACK_COLORS[tissue_key], linewidth=0.8, label=spec["label"])
    ax.set_xlabel("Position relative to rs12740374 (kb)")
    ax.set_ylabel("Synergy score (x1e3)")
    ax.legend(frameon=False, fontsize=7)
    _save_svg(fig, output)
