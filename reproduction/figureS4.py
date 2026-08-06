"""Clean-room reproduction for Figure S4 (A-D): sequencing depth and
frequency-weighted RNA predictions for Wang et al. 2018 CRISPR repair
outcomes in primary human hepatocytes.

Zero new AlphaGenome calls anywhere in this panel group:

- S4A is Wang read-count/composition statistics -- no AlphaGenome.
- S4B reuses Figure 2B's already-scored repair-product sequences. Verified
  before writing any scoring code: all 24 human-hepatocyte modeled indels
  are within Figure 2B's top-50 fused human+mouse selection, so calling
  `_score_sequences(..., panel="2B", ...)` here hits Figure 2B's existing
  cache and scores nothing new.
- S4C/S4D reuse Figure 2C's already-scored 11x11 systematic deletion grid.
  The design construction (`_build_deletion_grid`) uses the exact same
  CUT_LEFT/RS_POS/SEQ_LEN constants as Figure 2C, so calling
  `_score_sequences(..., panel="2C", ...)` here also hits an existing
  cache. The 30 prespecified designs (upstream 1-10 x downstream 0-2) are
  a subset of that grid; everything panel-specific (motif/protospacer/PAM
  overlap, post-edit junction reconstruction) is deterministic sequence-
  position arithmetic, not a new prediction.

S4A's total/unedited/other-indel read counts are the one exception to
"everything here is either fresh or reused": they come from a raw
sequencing-alignment summary Wang et al. do not redistribute in their
supplementary spreadsheet (the spreadsheet itself only lists the 24
already-called modeled indels, whose read counts already reproduce
`modeled_simple_indel_reads` exactly). Those three totals are therefore
frozen, documented constants, not fabricated -- matching this project's
established convention for external non-AlphaGenome numbers with no
traceable derivation script (`MANIFEST_NOTES.md`).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .common import Audit
from .figure1 import _save_svg
from .figure2 import (
    CHROM, CUT_LEFT, GENES, RS_MINOR, RS_POS, SEQ_LEN,
    _ag, _build_deletion_grid, _parse_wang, _replace_rs, _score_sequences, _wang_sequences,
    fetch_hg38, run_fig2b, run_fig2c, stage_wang,
)

MOTIF_START = 109_274_967
MOTIF_END = 109_274_976
MINOR_MOTIF = "GTTGCTCAAT"
PROTOSPACER_START = 109_274_963
PROTOSPACER_END = 109_274_982
PAM_START = 109_274_960
PAM_END = 109_274_962
TARGET_INDEX_S4B = 4  # matches Figure 2B's own default target-track index
TARGET_INDEX_S4CD = 3  # S4C/D use a different target-track index than Figure 2B/S4B

# Wang et al. 2018 report these primary-human-hepatocyte read totals in the
# main text; they are not present in, or derivable from, the redistributed
# supplementary spreadsheet (which only lists the 24 already-called modeled
# indels). See module docstring.
FROZEN_HUMAN_TOTAL_READS = 51_370
FROZEN_HUMAN_UNEDITED_READS = 49_932


def _overlap_length(start: int, end: int, feature_start: int, feature_end: int) -> int:
    """Overlap between two 1-based inclusive intervals."""
    return max(0, min(end, feature_end) - max(start, feature_start) + 1)


# --- S4A: Wang human-hepatocyte sequencing depth --------------------------

def run_figs4a(run_dir: Path, audit: Audit, *, wang_xls: Path | None = None) -> None:
    with audit.step("S4A: Wang human hepatocyte sequencing depth"):
        wang_path = stage_wang(run_dir, audit, wang_xls)
        human = _parse_wang(wang_path)
        human = human[human.experiment.eq("human_hepatocytes")]
        modeled_reads = int(human.reads.sum())
        if len(human) != 24:
            raise ValueError(f"Expected 24 human modeled indels, found {len(human)}")
        unedited = FROZEN_HUMAN_UNEDITED_READS
        total = FROZEN_HUMAN_TOTAL_READS
        total_indel_reads = total - unedited
        other = total_indel_reads - modeled_reads
        depth = pd.DataFrame([{
            "experiment": "Human hepatocytes",
            "total_reads": total,
            "unedited_reads": unedited,
            "modeled_simple_indel_reads": modeled_reads,
            "other_indel_reads": other,
            "editing_fraction_percent": 100.0 * total_indel_reads / total,
            "modeled_indel_coverage_percent": 100.0 * modeled_reads / total_indel_reads,
        }])
        out = run_dir / "derived/FigureS4A_sequencing_depth.tsv"
        out.parent.mkdir(parents=True, exist_ok=True)
        depth.to_csv(out, sep="\t", index=False)

        fig, ax = plt.subplots(figsize=(70 / 25.4, 45 / 25.4))
        row = depth.iloc[0]
        parts = [("Unedited", row.unedited_reads, "#D0D0D0"), ("Modeled simple indels", row.modeled_simple_indel_reads, "#4E79A7"), ("Other/complex indels", row.other_indel_reads, "#E15759")]
        left = 0.0
        for label, value, color in parts:
            ax.barh(0, value, left=left, color=color, label=label)
            left += value
        ax.set_yticks([])
        ax.set_xlabel("Reads")
        ax.legend(frameon=False, fontsize=6.5, loc="upper center", bbox_to_anchor=(0.5, -0.35), ncol=1)
        ax.text(0.02, 0.75, f"n = {row.total_reads:,}\n{row.editing_fraction_percent:.1f}% edited", transform=ax.transAxes, fontsize=6.5, va="top")
        _save_svg(fig, run_dir / "figures/FigureS4A.svg")


# --- S4B: frequency-weighted heatmap of all 24 human repair outcomes ------

def run_figs4b(run_dir: Path, audit: Audit, *, wang_xls: Path | None = None, fasta: Path | None = None) -> None:
    fasta_path = fasta or fetch_hg38(run_dir, audit)
    wang_path = stage_wang(run_dir, audit, wang_xls)
    indels, sequences = _wang_sequences(_parse_wang(wang_path), fasta_path)
    human = indels[indels.experiment.eq("human_hepatocytes")].sort_values("rank_in_sheet").reset_index(drop=True)
    if len(human) != 24:
        raise ValueError(f"Expected 24 human modeled indels, found {len(human)}")

    import hashlib

    import pysam

    from .figure2 import SequenceState

    # Reconstruct major/minor exactly as run_fig2b does, from the same fasta.
    start0 = (RS_POS - 1) - SEQ_LEN // 2
    rs_index = RS_POS - 1 - start0
    ref_seq = pysam.FastaFile(str(fasta_path)).fetch(CHROM, start0, start0 + SEQ_LEN).upper()
    minor_seq = _replace_rs(ref_seq, rs_index, RS_MINOR)
    major_hash = hashlib.sha256(ref_seq.encode("ascii")).hexdigest()
    minor_hash = hashlib.sha256(minor_seq.encode("ascii")).hexdigest()
    sequences.setdefault(major_hash, ref_seq)
    sequences.setdefault(minor_hash, minor_seq)

    states = [SequenceState("minor", minor_seq, minor_hash)]
    states.extend(SequenceState(f"edit_{digest[:16]}", sequences[digest], digest) for digest in human.sequence_sha256.astype(str))
    # Reusing Figure 2B's own cache (panel="2B"): every one of these
    # sequences is already scored there, so this makes no new API calls.
    rna = _score_sequences(states, run_dir, audit, "2B", batch_size=8, max_workers=4)

    target = rna[rna.target_index.astype(int).eq(TARGET_INDEX_S4B) & rna.gene_symbol.isin(GENES) & rna.tss_half_width.astype(int).eq(2000)].copy()
    baseline = target[target.sequence_sha256.eq(minor_hash)].set_index("gene_symbol").rna_mean_tss_pm2kb
    target["percent_change_edit_vs_minor"] = [
        100.0 * (value - baseline[gene]) / baseline[gene] for gene, value in zip(target.gene_symbol, target.rna_mean_tss_pm2kb, strict=False)
    ]

    by_hash = target[target.sequence_sha256.isin(human.sequence_sha256)].set_index(["gene_symbol", "sequence_sha256"]).percent_change_edit_vs_minor
    order = human.sequence_sha256.astype(str).tolist()
    raw = pd.DataFrame(index=GENES, columns=order, dtype=float)
    weighted = pd.DataFrame(index=GENES, columns=order, dtype=float)
    for gene in GENES:
        raw_pct = np.array([by_hash[(gene, digest)] for digest in order])
        raw.loc[gene] = raw_pct
        weighted.loc[gene] = human.p_i.to_numpy(dtype=float) * raw_pct

    def short_label(row: Any) -> str:
        sign = "+" if row.operation == "insertion" else "−"
        event_label = f"{sign}{row.event_sequence}" if int(row.event_length) <= 8 else f"{sign}{int(row.event_length)} bp"
        return f"{int(row.rank_in_sheet):02d} {event_label}"

    labels = [short_label(row) for row in human.itertuples(index=False)]
    raw.columns = labels
    weighted.columns = labels
    summary = pd.DataFrame({
        "gene": GENES,
        "frequency_weighted_mixture_percent": weighted.sum(axis=1).to_numpy(float),
        "n_human_products": weighted.shape[1],
        "modeled_frequency_sum": [float(human.p_i.sum())] * len(GENES),
        "target_index": [TARGET_INDEX_S4B] * len(GENES),
    })

    out = run_dir / "derived/FigureS4B_freqweighted_heatmap"
    out.mkdir(parents=True, exist_ok=True)
    weighted.reset_index(names="gene").to_csv(out / "matrix.csv", index=False)
    raw.reset_index(names="gene").to_csv(out / "unweighted_percent_matrix.csv", index=False)
    human.to_csv(out / "columns.csv", index=False)
    summary.to_csv(out / "summary.csv", index=False)
    _render_s4b(weighted, run_dir / "figures/FigureS4B.svg")


def _render_s4b(weighted: pd.DataFrame, output: Path) -> None:
    from matplotlib.colors import TwoSlopeNorm

    values = weighted.to_numpy(dtype=float)
    limit = float(np.ceil(np.nanmax(np.abs(values)) * 10.0) / 10.0) or 1.0
    norm = TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)
    fig, ax = plt.subplots(figsize=(145 / 25.4, 43 / 25.4))
    image = ax.pcolormesh(np.arange(weighted.shape[1] + 1), np.arange(weighted.shape[0] + 1), values, cmap="RdBu_r", norm=norm, edgecolors="white", linewidth=0.4)
    ax.set_aspect("equal")
    ax.set_ylim(weighted.shape[0], 0)
    ax.set_yticks(np.arange(weighted.shape[0]) + 0.5, weighted.index.tolist())
    ax.set_xticks(np.arange(weighted.shape[1]) + 0.5, weighted.columns.tolist(), rotation=90, fontsize=5.5)
    fig.colorbar(image, ax=ax, label="Weighted contribution (%)")
    _save_svg(fig, output)


# --- Shared: the 30 prespecified Figure 2C deletion designs ---------------

def _selected_grid_annotations(run_dir: Path, audit: Audit) -> pd.DataFrame:
    """The 30 designs (upstream 1-10 x downstream 0-2) from Figure 2C's grid,
    annotated with motif/protospacer/PAM overlap and post-edit junction
    reconstruction -- pure sequence-position arithmetic, no AlphaGenome."""
    fasta_path = fetch_hg38(run_dir, audit)
    designs, states, _minor_hash = _build_deletion_grid(fasta_path)
    selected = designs[
        designs.geometry.eq("full_xy_grid") & designs.operation.eq("deletion")
        & designs.upstream_bases.between(1, 10) & designs.downstream_bases.isin([0, 1, 2])
    ].copy()
    if len(selected) != 30:
        raise ValueError(f"Expected 30 selected designs, found {len(selected)}")
    selected["upstream_bases"] = selected["upstream_bases"].astype(int)
    selected["downstream_bases"] = selected["downstream_bases"].astype(int)
    selected = selected.sort_values(["downstream_bases", "upstream_bases"]).reset_index(drop=True)

    import pysam

    start0 = (RS_POS - 1) - SEQ_LEN // 2
    extra = 84
    extended = pysam.FastaFile(str(fasta_path)).fetch(CHROM, start0, start0 + SEQ_LEN + extra).upper()
    rs_index = RS_POS - 1 - start0
    cut_index = CUT_LEFT - start0
    minor = _replace_rs(extended, rs_index, RS_MINOR)

    rows = []
    for row in selected.itertuples(index=False):
        up, down = int(row.upstream_bases), int(row.downstream_bases)
        event_start = CUT_LEFT - up + 1
        event_end = CUT_LEFT + down
        deletion_start_idx = cut_index - up
        deletion_end_idx = cut_index + down
        deleted_sequence = minor[deletion_start_idx:deletion_end_idx]
        edited = minor[:deletion_start_idx] + minor[deletion_end_idx:]
        local_40bp = edited[cut_index - 20 : cut_index + 20]

        motif_overlap = _overlap_length(event_start, event_end, MOTIF_START, MOTIF_END)
        protospacer_overlap = _overlap_length(event_start, event_end, PROTOSPACER_START, PROTOSPACER_END)
        pam_overlap = _overlap_length(event_start, event_end, PAM_START, PAM_END)
        deletes_rs = event_start <= RS_POS <= event_end

        hits: list[int] = []
        search_from = 0
        while True:
            hit = local_40bp.find(MINOR_MOTIF, search_from)
            if hit < 0:
                break
            hits.append(hit)
            search_from = hit + 1
        # String index i of the (0-based) `minor` sequence is 1-based genomic
        # position start0+i+1; local_40bp's index 0 is string index
        # cut_index-20, so its genomic position is start0+(cut_index-20)+1,
        # not start0+(cut_index-20). Verified against a real local_40bp
        # window before trusting this (motif found at index 20, not 21).
        seq_local_start = start0 + (cut_index - 20) + 1
        junction_index = event_start - seq_local_start
        crossing_hits = [hit for hit in hits if hit < junction_index <= hit + len(MINOR_MOTIF) - 1]
        deletion_length = event_end - event_start + 1
        if event_start <= RS_POS <= event_end:
            rs_index_postedit = None
        elif RS_POS < event_start:
            rs_index_postedit = RS_POS - seq_local_start
        else:
            rs_index_postedit = RS_POS - seq_local_start - deletion_length

        if down == 0 and hits:
            status = "retained"
        elif crossing_hits:
            status = "junction_restored"
        elif rs_index_postedit is not None:
            status = "disrupted_rs_retained"
        else:
            status = "disrupted_rs_deleted"
        display_group = f"d{down}_{'intact' if down == 0 else ('restored' if status == 'junction_restored' else 'disrupted')}"

        rows.append({
            "design_id": row.design_id, "operation": row.operation, "length": row.length,
            "upstream_bases": up, "downstream_bases": down, "geometry": row.geometry,
            "event_start_hg38": event_start, "event_end_hg38_inclusive": event_end,
            "deletes_rs12740374": bool(deletes_rs), "deleted_sequence": deleted_sequence,
            "motif_overlap_bases": motif_overlap, "protospacer_overlap_bases": protospacer_overlap,
            "pam_overlap_bases": pam_overlap, "sequence_sha256": row.sequence_sha256,
            "local_sequence_40bp": local_40bp,
            "exact_minor_motif_present": bool(hits),
            "exact_minor_motif_hit_starts0": ";".join(str(v) for v in hits),
            "junction_crossing_motif_present": bool(crossing_hits),
            "junction_crossing_motif_hit_starts0": ";".join(str(v) for v in crossing_hits),
            "deletion_junction_index0": junction_index,
            "rs12740374_index_postedit0": rs_index_postedit,
            "postedit_motif_status": status, "display_group": display_group,
        })
    annotated = pd.DataFrame(rows)
    restored = {(int(r.upstream_bases), int(r.downstream_bases)) for r in annotated[annotated.postedit_motif_status.eq("junction_restored")].itertuples(index=False)}
    expected_restored = {(2, 1), (7, 1), (10, 1)}
    if restored != expected_restored:
        raise ValueError(f"Unexpected junction-restored designs: {sorted(restored)}")

    dna_state_by_hash = {s.sequence_sha256: s for s in states}
    scoring_states = [dna_state_by_hash[h] for h in annotated.sequence_sha256.unique() if h in dna_state_by_hash]
    minor_state = next(s for s in states if s.sequence_sha256 == _minor_hash)
    # Reusing Figure 2C's own cache (panel="2C"): every one of these
    # sequences is already scored there, so this makes no new API calls.
    rna = _score_sequences([minor_state, *scoring_states], run_dir, audit, "2C", batch_size=8, max_workers=4)
    target = rna[rna.target_index.astype(int).eq(TARGET_INDEX_S4CD) & rna.gene_symbol.isin(GENES) & rna.tss_half_width.astype(int).eq(2000)].copy()
    baseline = target[target.sequence_sha256.eq(_minor_hash)].set_index("gene_symbol").rna_mean_tss_pm2kb
    target["percent_change_vs_minor"] = [
        100.0 * (value - baseline[gene]) / baseline[gene] for gene, value in zip(target.gene_symbol, target.rna_mean_tss_pm2kb, strict=False)
    ]
    for gene in GENES:
        per_design = target[target.gene_symbol.eq(gene)].set_index("sequence_sha256").percent_change_vs_minor
        annotated[f"rna_change_percent_{gene}"] = annotated.sequence_sha256.map(per_design)
    return annotated


# --- S4C: junction reconstruction table ------------------------------------

def run_figs4c(run_dir: Path, audit: Audit) -> None:
    with audit.step("S4C: selected systematic deletion junction reconstruction"):
        annotated = _selected_grid_annotations(run_dir, audit)
        out = run_dir / "derived/FigureS4C_junction_reconstruction.tsv"
        out.parent.mkdir(parents=True, exist_ok=True)
        annotated.to_csv(out, sep="\t", index=False)
        _render_s4c(annotated, run_dir / "figures/FigureS4C.svg")


def _render_s4c(data: pd.DataFrame, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(120 / 25.4, 100 / 25.4))
    ordered = data.sort_values(["downstream_bases", "upstream_bases"]).reset_index(drop=True)
    status_color = {"retained": "#365f46", "junction_restored": "#c08a00", "disrupted_rs_retained": "#7a5a00", "disrupted_rs_deleted": "#a61b1b"}
    for i, row in ordered.iterrows():
        ax.barh(i, 1, color=status_color.get(row.postedit_motif_status, "#999999"))
    ax.set_yticks(range(len(ordered)), [f"u{r.upstream_bases:02d}/d{r.downstream_bases}" for r in ordered.itertuples()], fontsize=5.5)
    ax.set_xlim(0, 1)
    ax.set_xticks([])
    ax.set_title("Post-edit motif status")
    _save_svg(fig, output)


# --- S4D: downstream-stratified repair-class decomposition ----------------

def run_figs4d(run_dir: Path, audit: Audit) -> None:
    with audit.step("S4D: selected systematic deletion RNA effects by downstream extent"):
        annotated = _selected_grid_annotations(run_dir, audit)
        long_rows = []
        for gene in GENES:
            sub = annotated[["design_id", "upstream_bases", "downstream_bases", "display_group", f"rna_change_percent_{gene}"]].rename(
                columns={f"rna_change_percent_{gene}": "rna_change_percent"}
            )
            sub["gene_symbol"] = gene
            long_rows.append(sub)
        data = pd.concat(long_rows, ignore_index=True)
        if len(data) != 90:
            raise ValueError(f"Expected 90 rows (30 designs x 3 genes), found {len(data)}")
        summary = data.groupby(["gene_symbol", "display_group"], as_index=False).agg(
            n_designs=("design_id", "nunique"),
            median_rna_change_percent=("rna_change_percent", "median"),
            mean_rna_change_percent=("rna_change_percent", "mean"),
            min_rna_change_percent=("rna_change_percent", "min"),
            max_rna_change_percent=("rna_change_percent", "max"),
        )
        out = run_dir / "derived/FigureS4D_downstream_stratified"
        out.mkdir(parents=True, exist_ok=True)
        data.to_csv(out / "source_data.csv", index=False)
        summary.to_csv(out / "group_summary.csv", index=False)
        _render_s4d(data, run_dir / "figures/FigureS4D.svg")


def _render_s4d(data: pd.DataFrame, output: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(130 / 25.4, 55 / 25.4), sharey=True)
    for ax, gene in zip(axes, GENES, strict=True):
        sub = data[data.gene_symbol.eq(gene)]
        for downstream in (0, 1, 2):
            values = sub.loc[sub.downstream_bases.eq(downstream), "rna_change_percent"]
            ax.scatter(np.full(len(values), downstream), values, s=18, color="#3182bd", alpha=0.75)
            if len(values):
                ax.plot([downstream - 0.25, downstream + 0.25], [values.median()] * 2, color="#1f1f1f", linewidth=1.4)
        ax.axhline(0, color="#777777", lw=0.6)
        ax.set_xticks([0, 1, 2])
        ax.set_title(gene)
    axes[0].set_ylabel("RNA change vs minor (%)")
    fig.supxlabel("Bases deleted downstream of cut")
    _save_svg(fig, output)
