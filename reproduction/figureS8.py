"""Clean-room reproduction for Figure S8 (A-E): boundary robustness of the
315bp rs12740374 core module.

S8A (three-gene directional recovery) reuses Figure 3E's own design and
scoring cache directly (`_fig3e_build_design`/`_fig3e_retention`,
verified to use identical UPSTREAM_EXTENTS/DOWNSTREAM_EXTENTS/SEEDS
constants as the archive's own `run_single_arm_recovery.py`) -- zero new
AlphaGenome calls when `3E` was run in the same run directory first.

S8B (three-gene arm necessity, ALL_FOLDS vs FOLD_0) looks superficially
like Figure 3G's own component-necessity design (same four conditions:
downstream arm / upstream arm / C/EBP / both arms), and an earlier version
of this port reused Figure 3G's cache directly on that basis -- but a real
run showed a real, non-trivial mismatch (up to ~0.33 retention units, not
float noise) against the archive. Root cause, found by reading the
archive's actual generating script (`run_single_arm_scramble.py`, not
`run_component_necessity.py`/Figure 3G's source): the two use different
scramble constructions. Figure 3G's `_shuffle_groups` scrambles the
*complement* of a component directly from the native sequence. S8B's
`copy_template_positions` instead starts from a shared, once-per-seed
"armwise" scrambled template (`shared_armwise_scramble_template` --
already ported here as `_shared_armwise_scramble_template`, reused
byte-for-byte from Figure 3F, which independently needed the exact same
function) and copies scrambled positions *back into* a native background
for just the target arm interval -- the same "inside_scramble" logic
Figure 3F already implements in `_fig3f_nested_background`, just applied
to a single arm interval instead of a symmetric upstream/downstream
window. Ported here as `_arm_scramble_background`. Both ALL_FOLDS and
FOLD_0 are therefore genuinely fresh (~178 calls each, matching Figure
3G's own real cost for a similarly-sized design).

S8C-E (the exhaustive upstream x downstream boundary grid and its
seed-level derivatives) reuse the same `_shared_armwise_scramble_template`
plus Figure 3F's own `_fig3f_nested_background` and `_score_design(...,
panel="3F", ...)` cache directly -- any (upstream_bp, downstream_bp) pair
already inside Figure 3F's own 175-190 x 105-160 grid is a zero-cost hit;
only the pairs outside that rectangle are genuinely new. See
MANIFEST_NOTES.md for exact sizing.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .common import Audit, api_key
from .figure3 import (
    CHROM,
    GENES,
    LIVER,
    RS_POS,
    SEQ_LEN,
    _ag,
    _fetch_reference,
    _fig3e_build_design,
    _fig3e_retention,
    _score_design,
    _set_base,
    _shared_armwise_scramble_template,
    _summarize,
)

_S8B_LOCAL_HALF_BP = 1000
_S8B_SEEDS = list(range(1740374, 1740382))

# --- S8A: reuses Figure 3E's own design/cache -----------------------------

def run_figs8a(run_dir: Path, audit: Audit, fasta_path: Path, *, batch_size: int = 32, max_workers: int = 4) -> None:
    genome, _, _, _ = _ag()
    interval = genome.Interval(CHROM, RS_POS, RS_POS).resize(SEQ_LEN)
    ref_seq = _fetch_reference(fasta_path, interval)
    with audit.step("S8A: rebuild Figure 3E's directional-recovery design (reuses its own cache)"):
        sequences, design = _fig3e_build_design(ref_seq, int(interval.start))
    scored = _score_design(run_dir, audit, "3E", interval, sequences, design, batch_size=batch_size, max_workers=max_workers)
    with audit.step("S8A: aggregate per-gene arm-recovery summary"):
        retention = _fig3e_retention(scored)
        genes_only = retention[retention.gene.isin(GENES)].copy()
        summary_rows = []
        for (condition, arm, extent_bp, gene), group in genes_only.groupby(["condition", "arm", "extent_bp", "gene"]):
            deltas = group.delta_liver
            ret = group.retention_vs_native.dropna()
            summary_rows.append({
                "condition": condition, "arm": arm, "extent_bp": extent_bp, "gene": gene, "n_seeds": len(group),
                "mean_delta_liver": float(deltas.mean()), "sd_delta_liver": float(deltas.std()),
                "mean_retention": float(ret.mean()) if len(ret) else np.nan,
                "median_retention": float(ret.median()) if len(ret) else np.nan,
                "sd_retention": float(ret.std()) if len(ret) else np.nan,
                "min_retention": float(ret.min()) if len(ret) else np.nan,
                "max_retention": float(ret.max()) if len(ret) else np.nan,
                "positive_delta_fraction": float((deltas > 0).mean()),
            })
        summary = pd.DataFrame(summary_rows)
        out = run_dir / "derived/FigureS8_boundary_robustness"
        out.mkdir(parents=True, exist_ok=True)
        summary.to_csv(out / "S8A_summary.csv", index=False)
        _render_figs8a(summary, run_dir / "figures/FigureS8A.svg")


def _render_figs8a(summary: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(4.0, 2.6))
    colors = {"upstream_recovery": "#2478b5", "downstream_recovery": "#d62728"}
    for condition, color in colors.items():
        for gene in GENES:
            sub = summary[summary.condition.eq(condition) & summary.gene.eq(gene)].sort_values("extent_bp")
            ax.plot(sub.extent_bp, sub.mean_retention, color=color, alpha=0.4, lw=0.8)
    ax.axhline(1.0, color="#777777", ls=":", lw=0.85)
    ax.set_xlabel("Native sequence restored (bp)")
    ax.set_ylabel("Mean liver RNA retention")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)


# --- S8B: the archive's own "single-arm scramble" construction -----------
# (shared armwise-scramble template + arm-interval inside-scramble, NOT
# Figure 3G's component-necessity shuffle -- see module docstring).

_S8B_CONDITIONS: list[tuple[str, str, list[tuple[int, int]]]] = [
    ("upstream_only", "upstream", [(-179, -2)]),
    ("downstream_only", "downstream", [(9, 135)]),
    ("motif_only", "control", [(-1, -1), (1, 8)]),
    ("both_nonmotif_arms", "control", [(-179, -2), (9, 135)]),
    ("full_315", "control", [(-179, 135)]),
]


def _arm_scramble_background(ref_seq: str, scrambled_template: str, *, seq_start0: int, intervals: list[tuple[int, int]]) -> str:
    """The archive's `copy_template_positions`: native background, with
    the scrambled template's bases copied in over the given interval(s)
    (relative to rs12740374). The same "inside_scramble" direction as
    `_fig3f_nested_background`, generalized to arbitrary (not necessarily
    symmetric) intervals."""
    chars = list(ref_seq)
    rs_index = RS_POS - 1 - seq_start0
    for rel_start, rel_end in intervals:
        start_index, end_index = RS_POS + rel_start - 1 - seq_start0, RS_POS + rel_end - 1 - seq_start0
        for index in range(start_index, end_index + 1):
            if index != rs_index:
                chars[index] = scrambled_template[index]
    return "".join(chars)


def _figs8b_build_design(ref_seq: str, seq_start0: int) -> tuple[list[str], pd.DataFrame]:
    sequences: list[str] = []
    rows: list[dict[str, object]] = []
    for seed in _S8B_SEEDS:
        template = _shared_armwise_scramble_template(ref_seq, seq_start0=seq_start0, local_half_bp=_S8B_LOCAL_HALF_BP, scramble_seed=seed)
        for condition, arm, intervals in [("native", "none", [])] + _S8B_CONDITIONS:
            background = ref_seq if condition == "native" else _arm_scramble_background(ref_seq, template, seq_start0=seq_start0, intervals=intervals)
            for allele, base_value in (("REF", "G"), ("ALT", "T")):
                sequence = _set_base(background, seq_start0, RS_POS, base_value)
                sequences.append(sequence)
                rows.append({
                    "sequence_index": len(sequences) - 1, "scramble_seed": seed, "condition": condition, "arm": arm,
                    "allele": allele, "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                })
    return sequences, pd.DataFrame(rows)


def _figs8b_retention(scored: pd.DataFrame) -> pd.DataFrame:
    keys = ["scramble_seed", "condition", "arm", "gene"]
    ref = scored[scored.allele.eq("REF")].copy()
    alt = scored[scored.allele.eq("ALT")].copy()
    paired = ref.merge(alt, on=keys, suffixes=("_ref", "_alt"), validate="one_to_one")
    paired["delta_liver"] = paired.liver_rna_signal_alt - paired.liver_rna_signal_ref
    native = paired[paired.condition.eq("native")][["scramble_seed", "gene", "delta_liver"]].rename(columns={"delta_liver": "native_delta_liver"})
    paired = paired.merge(native, on=["scramble_seed", "gene"], validate="many_to_one")
    denominator = paired.native_delta_liver.abs().where(paired.native_delta_liver.abs() > 1e-12, np.nan)
    paired["retention"] = paired.delta_liver.abs() / denominator
    return paired


def _score_figs8b_design(run_dir: Path, audit: Audit, panel: str, model_version_name: str, interval: Any, sequences: list[str], design: pd.DataFrame, *, batch_size: int, max_workers: int) -> pd.DataFrame:
    genome, dna_client, dna_model, _ = _ag()
    cache = run_dir / "predictions" / f"Figure{panel}_sequence_cache"
    cache.mkdir(parents=True, exist_ok=True)
    hash_to_seq = dict(zip(design.sequence_sha256, [sequences[int(i)] for i in design.sequence_index]))
    values: dict[str, dict[str, float]] = {}
    for digest in hash_to_seq:
        path = cache / f"{digest}.tsv"
        if path.exists():
            values[digest] = pd.read_csv(path, sep="\t").set_index("gene").liver_rna_signal.to_dict()
    missing = [digest for digest in hash_to_seq if digest not in values]
    client = None
    for start in range(0, len(missing), batch_size):
        batch = missing[start : start + batch_size]
        if client is None:
            client = dna_client.create(api_key(), model_version=getattr(dna_model.ModelVersion, model_version_name), timeout=300)
        with audit.step(f"{panel}: score sequences {start + 1}-{start + len(batch)} ({model_version_name})"):
            outputs = client.predict_sequences(
                sequences=[hash_to_seq[digest] for digest in batch],
                requested_outputs={dna_client.OutputType.RNA_SEQ}, ontology_terms=[LIVER],
                intervals=[interval] * len(batch), progress_bar=False, max_workers=max_workers,
            )
            for digest, output in zip(batch, outputs, strict=True):
                gene_values = _summarize(output, interval)
                pd.DataFrame([{"gene": g, "liver_rna_signal": gene_values[g]} for g in GENES]).to_csv(cache / f"{digest}.tsv", sep="\t", index=False)
                values[digest] = gene_values
            audit.add_api_calls(panel, len(batch))
            audit.add_api_requests(panel, 1)
    rows: list[dict[str, object]] = []
    for row in design.itertuples(index=False):
        gene_values = values[row.sequence_sha256]
        for gene in GENES:
            rows.append({**row._asdict(), "gene": gene, "liver_rna_signal": gene_values[gene]})
    return pd.DataFrame(rows)


def _figs8b_summary(retention: pd.DataFrame) -> pd.DataFrame:
    """Same full per-condition/arm/gene summary schema as S8A (n_seeds,
    mean/sd delta, mean/median/sd/min/max retention, positive_delta_fraction)."""
    genes_only = retention[retention.gene.isin(GENES)].copy()
    rows = []
    conditions = [("native", "none")] + [(c, a) for c, a, _ in _S8B_CONDITIONS]
    for condition, arm in conditions:
        sub = genes_only[genes_only.condition.eq(condition)]
        for gene, group in sub.groupby("gene"):
            deltas = group.delta_liver
            ret = group.retention.dropna()
            rows.append({
                "condition": condition, "arm": arm, "extent_bp": np.nan, "gene": gene, "n_seeds": len(group),
                "mean_delta_liver": float(deltas.mean()), "sd_delta_liver": float(deltas.std()),
                "mean_retention": float(ret.mean()) if len(ret) else np.nan,
                "median_retention": float(ret.median()) if len(ret) else np.nan,
                "sd_retention": float(ret.std()) if len(ret) else np.nan,
                "min_retention": float(ret.min()) if len(ret) else np.nan,
                "max_retention": float(ret.max()) if len(ret) else np.nan,
                "positive_delta_fraction": float((deltas > 0).mean()),
            })
    return pd.DataFrame(rows)


def run_figs8b(run_dir: Path, audit: Audit, fasta_path: Path, *, batch_size: int = 32, max_workers: int = 4) -> None:
    genome, _, _, _ = _ag()
    interval = genome.Interval(CHROM, RS_POS, RS_POS).resize(SEQ_LEN)
    ref_seq = _fetch_reference(fasta_path, interval)
    with audit.step("S8B: build the single-arm scramble design (shared armwise template)"):
        sequences, design = _figs8b_build_design(ref_seq, int(interval.start))

    scored_all_folds = _score_figs8b_design(run_dir, audit, "S8B_all_folds", "ALL_FOLDS", interval, sequences, design, batch_size=batch_size, max_workers=max_workers)
    scored_fold0 = _score_figs8b_design(run_dir, audit, "S8B_fold0", "FOLD_0", interval, sequences, design, batch_size=batch_size, max_workers=max_workers)

    with audit.step("S8B: aggregate arm-necessity summary (ALL_FOLDS vs FOLD_0)"):
        summary_all_folds = _figs8b_summary(_figs8b_retention(scored_all_folds))
        summary_fold0 = _figs8b_summary(_figs8b_retention(scored_fold0))
        out = run_dir / "derived/FigureS8_boundary_robustness"
        out.mkdir(parents=True, exist_ok=True)
        summary_all_folds.to_csv(out / "S8B_arm_audit_all_folds_summary.csv", index=False)
        summary_fold0.to_csv(out / "S8B_arm_audit_fold0_summary.csv", index=False)
        _render_figs8b(summary_all_folds, summary_fold0, run_dir / "figures/FigureS8B.svg")


def _render_figs8b(summary_all_folds: pd.DataFrame, summary_fold0: pd.DataFrame, path: Path) -> None:
    display = {"downstream_only": "Downstream\narm", "upstream_only": "Upstream\narm", "motif_only": "C/EBP", "both_nonmotif_arms": "Both\narms"}
    order = list(display)
    x = np.arange(len(order))
    fig, axes = plt.subplots(1, len(GENES), figsize=(6.5, 2.4), sharey=True)
    for ax, gene in zip(axes, GENES, strict=True):
        for model, summary, color, marker, shift in [
            ("ALL_FOLDS", summary_all_folds, "#2478b5", "o", -0.09), ("FOLD_0", summary_fold0, "#e68613", "D", 0.09),
        ]:
            values = summary[summary.gene.eq(gene)].set_index("condition")
            means = [values.loc[c, "mean_retention"] if c in values.index else np.nan for c in order]
            sems = [values.loc[c, "sd_retention"] / np.sqrt(values.loc[c, "n_seeds"]) if c in values.index else np.nan for c in order]
            ax.errorbar(x + shift, means, yerr=sems, color=color, marker=marker, ms=4, lw=0, elinewidth=0.9, capsize=1.8, label=model)
        ax.axhline(1.0, color="#777777", ls=":", lw=0.85)
        ax.set_xticks(x, [display[c] for c in order])
        ax.set_title(gene, fontsize=8)
    axes[0].set_ylabel("Liver RNA retention")
    axes[1].legend(frameon=False, fontsize=6.5, ncol=2, loc="upper center")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)
