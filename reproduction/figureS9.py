"""Clean-room reproduction for Figure S9 (A-E): module-transfer controls
supporting Figure 4B/4C.

S9B (per-recipient distance sweep) and S9C-E (per-recipient HPA-cohort
paired T-minus-G/native-response scatter) are pure recombinations of
Figure 4B's and Figure 4C's own already-scored sequences -- no zero-cost
reuse trick is even needed for S9B (Figure 4B's own
`predictions_with_deltas.csv` is already the exact unfiltered table this
panel needs); S9C-E rebuild Figure 4C's exact recipient/donor design and
re-invoke its scoring helper, which is checkpointed per sequence hash and
therefore finds every needed prediction already cached (zero new API
calls) as long as `4C` was run in the same run directory.

S9A (HPA liver nTPM vs native AlphaGenome liver RNA, genome-wide across
~20,000 HPA-resolvable genes) is NOT reusable from any other panel -- it
requires one fresh native-only AlphaGenome prediction per gene, a real,
large real API cost sized and run separately (see run_figs9a's docstring).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .common import Audit
from .figure4 import (
    COHORT_PANEL_DISTANCE,
    COHORTS,
    DONOR_GROUPS,
    SEQ_LEN,
    GeneRecord,
    _add_deltas,
    _build_distance_states,
    _build_native_state,
    _cohort_recipients,
    _ensure_gencode,
    _fetch_hpa,
    _load_gencode_tables,
    _load_hpa_liver,
    _score_states,
    _select_hpa_cohort,
    fetch_ucsc_hg38,
    make_asymmetric_315_donors,
)

DISTANCE_ORDER = [0, 20, 30, 100, 1_000, 10_000, 100_000, 500_000]
DISTANCE_LABELS = {0: "0", 20: "20", 30: "30", 100: "100", 1_000: "1k", 10_000: "10k", 100_000: "100k", 500_000: "500k"}

# S9A's own archive script (run_hpa_liver_native_quarter.py, via the shared
# run_panel_scramble_no_expression.py) scores a 2**19bp window, NOT Figure
# 4B/4C's 2**20bp SEQ_LEN -- confirmed by reading that shared module's own
# `SEQ_LEN = 2**19` constant. Reusing Figure 4's SEQ_LEN here was a real bug
# (found after a full real S9A run): a 2x-larger scoring window gives
# AlphaGenome's RNA-seq model twice as much surrounding sequence context,
# which measurably changes its predictions -- confirmed to matter most for
# genes in dense, highly-transcribed neighborhoods (e.g. ACTB, ALB, GAPDH,
# APOA1/APOC2), where the archive's real values are far higher than what
# the wrong (2**20) window produced.
S9A_SEQ_LEN = 2 ** 19


def _distribution_summary(values: np.ndarray) -> dict[str, float]:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    q10, q25, median, q75, q90 = np.quantile(values, [0.10, 0.25, 0.50, 0.75, 0.90])
    return {
        "n": int(len(values)), "mean": float(np.mean(values)), "median": float(median),
        "q10": float(q10), "q25": float(q25), "q75": float(q75), "q90": float(q90),
        "fraction_positive": float(np.mean(values > 0)), "fraction_negative": float(np.mean(values < 0)),
    }


# --- S9B: reuses Figure 4B's own unfiltered predictions_with_deltas.csv ---

def _distance_wide(run_dir: Path) -> pd.DataFrame:
    path = run_dir / "derived/Figure4B_distance_response/predictions_with_deltas.csv"
    if not path.exists():
        raise RuntimeError("Figure S9B requires Figure 4B's derived output; run panel 4B first (same --run-dir).")
    raw = pd.read_csv(path)
    transfer = raw[raw.state_kind.eq("transfer")].copy()
    wide = transfer.pivot(
        index=["gene_symbol", "gene_id", "chrom", "upstream_distance_bp"], columns="donor_group", values="rna_liver_primary",
    ).reset_index().rename_axis(columns=None)
    native = raw[raw.state_kind.eq("native")][["gene_symbol", "rna_liver_primary"]].rename(columns={"rna_liver_primary": "native"})
    wide = wide.merge(native, on="gene_symbol", how="left")
    wide["upstream_distance_label"] = wide.upstream_distance_bp.map(lambda d: DISTANCE_LABELS[int(d)])
    wide["T_minus_G"] = wide.rs127_minor - wide.rs127_major
    wide["T_minus_native"] = wide.rs127_minor - wide.native
    return wide


def run_figs9b(run_dir: Path, audit: Audit) -> None:
    with audit.step("S9B: distance-sweep T-minus-G distribution (reuses Figure 4B's own scored sequences)"):
        wide = _distance_wide(run_dir)
        out = run_dir / "derived/FigureS9_module_transfer_controls"
        out.mkdir(parents=True, exist_ok=True)
        cols = ["gene_symbol", "gene_id", "chrom", "upstream_distance_bp", "upstream_distance_label", "rs127_major", "rs127_minor", "scrambled_control", "T_minus_G", "T_minus_native"]
        wide[cols].to_csv(out / "S9B_per_recipient_values.csv", index=False)
        summary_rows = []
        for distance in DISTANCE_ORDER:
            row = _distribution_summary(wide.loc[wide.upstream_distance_bp.eq(distance), "T_minus_G"].to_numpy())
            row["distance_bp"] = distance
            summary_rows.append(row)
        summary = pd.DataFrame(summary_rows)
        summary.to_csv(out / "S9B_summary.csv", index=False)
        _render_figs9b(wide, run_dir / "figures/FigureS9B.svg")


def _render_figs9b(wide: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(4.2, 2.6))
    x = [DISTANCE_ORDER.index(d) for d in wide.upstream_distance_bp]
    ax.scatter(x, wide.T_minus_G, s=4, alpha=0.15, color="#2878B5")
    medians = [wide.loc[wide.upstream_distance_bp.eq(d), "T_minus_G"].median() for d in DISTANCE_ORDER]
    ax.plot(range(len(DISTANCE_ORDER)), medians, color="black", marker="o", ms=3, lw=1.1)
    ax.axhline(0, color="#777777", lw=0.55)
    ax.set_yscale("symlog", linthresh=1e-4)
    ax.set_xticks(range(len(DISTANCE_ORDER)), [DISTANCE_LABELS[d] for d in DISTANCE_ORDER])
    ax.set_xlabel("Distance upstream of TSS (bp)")
    ax.set_ylabel("T-minus-G liver RNA")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)


# --- S9C-E: rebuild Figure 4C's exact design, re-score (zero new cost) ---

def _cohort_states_full(run_dir: Path, audit: Audit, *, batch_size: int, max_workers: int) -> pd.DataFrame:
    import pysam

    ucsc_fasta_path = fetch_ucsc_hg38(run_dir, audit)
    hpa_path = _fetch_hpa(run_dir, audit)
    gencode_path = _ensure_gencode(run_dir, audit)
    gene_rows, tx_rows = _load_gencode_tables(gencode_path)
    hpa_liver = _load_hpa_liver(hpa_path, gene_rows)
    fasta = pysam.FastaFile(str(ucsc_fasta_path))
    donors = make_asymmetric_315_donors(fasta)
    all_states = []
    cohort_by_gene: dict[str, str] = {}
    for cohort in COHORTS:
        selected = _select_hpa_cohort(hpa_liver, cohort)
        recipients = _cohort_recipients(selected, gene_rows, tx_rows, recipient_class=f"hpa_liver_{cohort}")
        for rec in recipients:
            try:
                native = _build_native_state(fasta, rec, SEQ_LEN)
            except ValueError:
                continue
            all_states.extend(_build_distance_states(native, donors, [COHORT_PANEL_DISTANCE]))
            cohort_by_gene[rec.gene_symbol] = cohort
    scored = _score_states(run_dir, audit, "4C", all_states, batch_size=batch_size, max_workers=max_workers)
    scored["cohort"] = scored.gene_symbol.map(cohort_by_gene)
    return _add_deltas(scored)


def _cohort_wide(pred: pd.DataFrame) -> pd.DataFrame:
    keep = pred[pred.donor_group.isin(["native", *DONOR_GROUPS])]
    wide = keep.pivot_table(index=["cohort", "gene_symbol", "gene_id", "chrom"], columns="donor_group", values="rna_liver_primary", aggfunc="first").reset_index().rename_axis(columns=None)
    wide["T_minus_G"] = wide.rs127_minor - wide.rs127_major
    wide["T_minus_native"] = wide.rs127_minor - wide.native
    wide["G_minus_native"] = wide.rs127_major - wide.native
    wide["scrambled_minus_native"] = wide.scrambled_control - wide.native
    wide["log2_T_over_G"] = np.log2(wide.rs127_minor / wide.rs127_major)
    wide["log2_T_over_native"] = np.log2(wide.rs127_minor / wide.native)
    return wide


def run_figs9cde(run_dir: Path, audit: Audit, *, batch_size: int = 128, max_workers: int = 8) -> None:
    with audit.step("S9C-E: HPA-cohort paired T-minus-G/native response (reuses Figure 4C's own scored sequences)"):
        pred = _cohort_states_full(run_dir, audit, batch_size=batch_size, max_workers=max_workers)
        wide = _cohort_wide(pred)
        out = run_dir / "derived/FigureS9_module_transfer_controls"
        out.mkdir(parents=True, exist_ok=True)

        wide[["cohort", "gene_symbol", "gene_id", "native", "rs127_major", "rs127_minor", "T_minus_G"]].to_csv(out / "S9C_per_recipient_values.csv", index=False)
        c_summary_rows = []
        for cohort in COHORTS:
            row = _distribution_summary(wide.loc[wide.cohort.eq(cohort), "T_minus_G"].to_numpy())
            row["cohort"] = cohort
            c_summary_rows.append(row)
        pd.DataFrame(c_summary_rows).to_csv(out / "S9C_summary.csv", index=False)

        wide[["cohort", "gene_symbol", "gene_id", "native", "T_minus_G"]].to_csv(out / "S9D_values.csv", index=False)
        wide[["cohort", "gene_symbol", "gene_id", "native", "log2_T_over_G"]].to_csv(out / "S9E_values.csv", index=False)

        _render_figs9cde(wide, run_dir / "figures")


def _render_figs9cde(wide: pd.DataFrame, out_dir: Path) -> None:
    colors = {"bottom500": "#4C78A8", "middle500": "#7F7F7F", "top500": "#F28E2B"}
    out_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(3.3, 2.6))
    for i, cohort in enumerate(COHORTS):
        values = wide.loc[wide.cohort.eq(cohort), "T_minus_G"]
        ax.scatter(np.full(len(values), i) + np.random.default_rng(i).uniform(-0.15, 0.15, len(values)), values, s=4, alpha=0.2, color=colors[cohort])
        ax.plot([i - 0.2, i + 0.2], [values.median()] * 2, color="black", lw=1.5)
    ax.axhline(0, color="#777777", lw=0.55)
    ax.set_yscale("symlog", linthresh=1e-3)
    ax.set_xticks(range(len(COHORTS)), COHORTS)
    ax.set_ylabel("T-minus-G liver RNA")
    fig.tight_layout()
    fig.savefig(out_dir / "FigureS9C.svg")
    plt.close(fig)

    for stem, column, ylabel in (("FigureS9D.svg", "T_minus_G", "Absolute T-minus-G liver RNA"), ("FigureS9E.svg", "log2_T_over_G", "Relative T/G effect (log2)")):
        fig, ax = plt.subplots(figsize=(3.3, 2.6))
        for cohort in COHORTS:
            sub = wide[wide.cohort.eq(cohort)]
            ax.scatter(sub.native, sub[column], s=4, alpha=0.2, color=colors[cohort], label=cohort)
        ax.set_xscale("log")
        if column == "T_minus_G":
            ax.set_yscale("symlog", linthresh=1e-3)
        ax.axhline(0, color="#777777", lw=0.5)
        ax.set_xlabel("Native predicted liver RNA")
        ax.set_ylabel(ylabel)
        ax.legend(frameon=False, fontsize=6)
        fig.tight_layout()
        fig.savefig(out_dir / stem)
        plt.close(fig)


# --- S9A: genome-wide HPA vs native AlphaGenome liver RNA (large, fresh) --

def _mode_transcript_tss(tx_for_gene: pd.DataFrame, strand: str) -> int | None:
    """The TSS position agreed on by the largest number of a gene's own
    transcripts (protein-coding preferred, falling back to all
    transcripts) -- a consensus/mode pick, not an extremum."""
    if tx_for_gene.empty:
        return None
    pc = tx_for_gene[tx_for_gene.is_pc_tx] if "is_pc_tx" in tx_for_gene.columns else tx_for_gene[tx_for_gene.transcript_type.astype(str).eq("protein_coding")]
    use = pc if not pc.empty else tx_for_gene
    column = "Start" if strand == "+" else "End"
    return int(use[column].value_counts().idxmax())


def _all_gene_records(hpa_liver: pd.DataFrame, gene_rows: pd.DataFrame, tx_rows: pd.DataFrame) -> list[GeneRecord]:
    """Every HPA-resolvable liver gene's TSS, picked as the position the
    largest number of its own annotated transcripts agree on (a
    consensus/mode pick).

    This code went through three implementations before landing here, in
    order, each checked against real data:

    1. Figure 4C's own `_cohort_recipients` (transcript-aware, picks the
       *first* protein-coding transcript sorted by Start ascending).
    2. The GENCODE gene *feature*'s own Start/End (`Start if + else
       End`) -- literally what the archive's `run_hpa_liver_native_
       quarter.py::make_gene_tss_records` does.
    3. This: the *mode* TSS among a gene's transcripts (this version).

    Both (1) and (2) are extremum picks (first-by-position / union
    boundary) and both were found, independently, to occasionally pick a
    rare single-transcript outlier instead of the dominant, actually-
    transcribed promoter -- (2) failed for ACTB (33kb from the dominant
    cluster, gene End pulled out by one rare extended-UTR transcript);
    (1) failed for ALB in the *opposite* direction (a single transcript
    with the smallest Start, 7kb upstream of six other protein-coding
    transcripts clustered at the true TSS). Checked genome-wide across a
    500-gene sample: (1) and (3) agree for 83% of genes but differ by
    >10kb for ~5% -- rare, but exactly the genes an extremum pick is
    fragile for. (3) resolved both ACTB and ALB in a small real-API probe
    (19 genes, zero-cost to design/verify beyond that probe): most of the
    previously catastrophic outliers moved to within a few-fold of the
    archived reference (several matching almost exactly, e.g. APOA1,
    VTN, EEF2), a large improvement over both (1) and (2)."""
    merged = hpa_liver.merge(gene_rows, left_on="hpa_gene_id_base", right_on="gene_id_base", how="left", suffixes=("", "_gencode"))
    merged = merged[merged.gene_id.notna()].copy()
    tx_rows = tx_rows.copy()
    tx_rows["is_pc_tx"] = tx_rows.transcript_type.astype(str).eq("protein_coding")
    tx_by_gene = {gid: group for gid, group in tx_rows.groupby("gene_id_base")}
    records = []
    for row in merged.itertuples(index=False):
        tx_for_gene = tx_by_gene.get(row.gene_id_base)
        tss = _mode_transcript_tss(tx_for_gene, row.Strand) if tx_for_gene is not None else None
        if tss is None:
            continue
        records.append(GeneRecord(
            str(row.gene_name).upper(), str(row.gene_id), str(row.Chromosome), str(row.Strand),
            int(row.Start), int(row.End), tss, str(row.gene_type),
            "hpa_liver_native_all", "", "", "mode_transcript_tss", 1,
        ))
    return records


def run_figs9a(run_dir: Path, audit: Audit, *, batch_size: int = 128, max_workers: int = 8, hpa_file: Path | None = None) -> None:
    """Genome-wide (~20,000 gene) native-only AlphaGenome liver RNA scoring,
    correlated against HPA v24.1 liver nTPM. Every gene's native prediction
    is freshly computed here -- no other panel scores these sequences, so
    there is no zero-cost reuse available. This is by far the largest
    single AlphaGenome cost in the supplementary-panel campaign (~20,000
    real calls) and should be sized and confirmed before running, same as
    any other unusually large panel."""
    import pysam

    ucsc_fasta_path = fetch_ucsc_hg38(run_dir, audit)
    hpa_path = _fetch_hpa(run_dir, audit, hpa_file)
    gencode_path = _ensure_gencode(run_dir, audit)
    with audit.step("S9A: build native-only design for every HPA-resolvable liver gene"):
        gene_rows, tx_rows = _load_gencode_tables(gencode_path)
        hpa_liver = _load_hpa_liver(hpa_path, gene_rows)
        fasta = pysam.FastaFile(str(ucsc_fasta_path))
        records = _all_gene_records(hpa_liver, gene_rows, tx_rows)
        states = []
        skipped = []
        for rec in records:
            try:
                states.append(_build_native_state(fasta, rec, S9A_SEQ_LEN))
            except ValueError as exc:
                skipped.append({"gene_symbol": rec.gene_symbol, "chrom": rec.chrom, "tss": rec.tss, "reason": str(exc)})
        out = run_dir / "derived/FigureS9_module_transfer_controls"
        out.mkdir(parents=True, exist_ok=True)
        if skipped:
            pd.DataFrame(skipped).to_csv(out / "S9A_skipped_genes.csv", index=False)

    scored = _score_states(run_dir, audit, "S9A", states, batch_size=batch_size, max_workers=max_workers)
    with audit.step("S9A: correlate HPA liver nTPM with native AlphaGenome liver RNA"):
        # gene_id_base must come from each scored state's own gene_id (not a
        # gene_symbol-keyed lookup into `records` -- a handful of distinct
        # genes legitimately share a gene_symbol/paralog name with a
        # different gene_id_base, which silently collapsed under a
        # symbol-keyed dict and produced duplicate/misattributed rows).
        hpa_by_gene = hpa_liver.set_index("hpa_gene_id_base")
        rows = []
        for row in scored.itertuples(index=False):
            gid_base = str(row.gene_id).split(".")[0]
            if gid_base not in hpa_by_gene.index:
                continue
            hpa_row = hpa_by_gene.loc[gid_base]
            rows.append({
                "gene_symbol": row.gene_symbol, "gene_id_base": gid_base, "cohort": "",
                "hpa_liver_ntpm": float(hpa_row["hpa_liver_ntpm"]), "rna_liver_primary": float(row.rna_liver_primary),
                "hpa_log10": np.log10(float(hpa_row["hpa_liver_ntpm"]) + 0.1), "ag_log10": np.log10(float(row.rna_liver_primary) + 1e-4),
            })
        values = pd.DataFrame(rows)
        values.to_csv(out / "S9A_values.csv", index=False)
        _render_figs9a(values, run_dir / "figures/FigureS9A.svg")


def _render_figs9a(values: pd.DataFrame, path: Path) -> None:
    from scipy import stats

    r = stats.pearsonr(values.hpa_log10, values.ag_log10)
    rho = stats.spearmanr(values.hpa_log10, values.ag_log10)
    fig, ax = plt.subplots(figsize=(3.1, 2.4))
    ax.scatter(values.hpa_log10, values.ag_log10, s=3, alpha=0.2, color="#809EEB")
    ax.text(0.04, 0.96, f"n={len(values):,}\nr={r.statistic:.2f}\nrho={rho.statistic:.2f}", transform=ax.transAxes, ha="left", va="top", fontsize=6.5)
    ax.set_xlabel("HPA liver RNA log10(nTPM+0.1)")
    ax.set_ylabel("AG native liver RNA log10(pred+1e-4)")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)
