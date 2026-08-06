"""Clean-room reproduction for Figure S6 (A-C): locus-matched held-out
replication of the Figure 2F multi-element Kircher accessibility benchmark,
restricted to the five non-SORT1 elements (SORT1 itself is covered by
Figure S5).

The ALL_FOLDS pass points `_score_multielement` at Figure 2F's own
`Figure2F_kircher_multielement` cache directory -- `load_multielement`
produces the identical experimental table and scan range Figure 2F already
scores, so if `2F` was run earlier in the *same* run directory, every chunk
is already on disk and the ALL_FOLDS pass costs zero new API calls (run
e.g. `--panels 2F,S6A,S6B,S6C` to get this). A standalone `--panels
S6A,S6B,S6C` run (no prior 2F in that run directory) finds no cache and
scores ALL_FOLDS fresh instead -- still a fully valid, freshly-computed
result, just not free. Either way, each element's own genomic-fold-matched
held-out pass is always genuinely new (Figure 2F never touches it). The
matched models are FOLD_0 for F9, FOLD_1 for FOXE1 and MYC, FOLD_2 for
PKLR, and FOLD_3 for LDLR (manuscript methods).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy import stats

from .common import Audit, api_key
from .figure2 import ELEMENTS, _ag, _score_multielement, load_multielement

OTHER_ELEMENTS = ("F9", "FOXE1", "LDLR", "MYC_rs11986220", "PKLR")
HELDOUT_FOLD = {"F9": "FOLD_0", "FOXE1": "FOLD_1", "LDLR": "FOLD_3", "MYC_rs11986220": "FOLD_1", "PKLR": "FOLD_2"}
DISPLAY = {"F9": "F9", "FOXE1": "FOXE1", "LDLR": "LDLR", "MYC_rs11986220": "MYC", "PKLR": "PKLR"}
MODEL_COLORS = {"ALL_FOLDS": "#377EB8", "matched held-out": "#E68632"}

SOURCE_COLUMNS = [
    "Chrom", "Pos", "Ref", "Alt", "mpra_effect", "mpra_sd", "mpra_n_experiments",
    "mpra_min_barcodes", "mpra_min_p", "element", "ag_accessibility_effect",
    "n_matched_tracks", "ontology", "modality", "model", "model_version",
]


def run_figs6(run_dir: Path, audit: Audit, kircher_path: Path) -> None:
    _, dna_client, dna_model, _ = _ag()
    experimental = load_multielement(kircher_path)
    parts = []
    for element in OTHER_ELEMENTS:
        table = experimental[element]

        # ALL_FOLDS: identical data/range to Figure 2F -- reuses its cache.
        with audit.step(f"S6: score {element} accessibility (ALL_FOLDS, reused from 2F)"):
            client = dna_client.create(api_key(), model_version=dna_model.ModelVersion.ALL_FOLDS, timeout=300)
            for modality in ELEMENTS[element][2]:
                scored = _score_multielement(
                    client, element, table, modality, run_dir, audit,
                    panel="S6_all_folds", cache_dirname="Figure2F_kircher_multielement",
                )
                merged = table.merge(scored, on=["Chrom", "Pos", "Ref", "Alt"], how="inner")
                merged["model"] = "ALL_FOLDS"
                merged["model_version"] = ""
                parts.append(merged)

        # Matched held-out fold: genuinely new, Figure 2F never scores this.
        fold = HELDOUT_FOLD[element]
        with audit.step(f"S6: score {element} accessibility (matched held-out {fold})"):
            client = dna_client.create(api_key(), model_version=getattr(dna_model.ModelVersion, fold), timeout=300)
            for modality in ELEMENTS[element][2]:
                scored = _score_multielement(
                    client, element, table, modality, run_dir, audit,
                    panel="S6_matched_heldout", cache_dirname=f"FigureS6_{element}_{fold}_heldout",
                )
                merged = table.merge(scored, on=["Chrom", "Pos", "Ref", "Alt"], how="inner")
                merged["model"] = "matched held-out"
                merged["model_version"] = fold
                parts.append(merged)

    merged_all = pd.concat(parts, ignore_index=True)[SOURCE_COLUMNS]
    counts = merged_all.groupby(["element", "modality", "model"]).size().unstack("model")
    if not (counts["ALL_FOLDS"] == counts["matched held-out"]).all():
        raise ValueError("ALL_FOLDS and held-out comparisons do not contain identical edits")

    statistics = _element_statistics(merged_all)

    out_dir = run_dir / "derived" / "FigureS6_kircher_other_loci"
    out_dir.mkdir(parents=True, exist_ok=True)
    dnase = merged_all[merged_all.modality.eq("DNASE")]
    atac = merged_all[merged_all.modality.eq("ATAC")]
    dnase.to_csv(out_dir / "FigureS6A_DNase_ALL_FOLDS_vs_matched_heldout_source.tsv", sep="\t", index=False)
    atac.to_csv(out_dir / "FigureS6B_ATAC_ALL_FOLDS_vs_matched_heldout_source.tsv", sep="\t", index=False)
    statistics.to_csv(out_dir / "FigureS6C_accessibility_correlation_model_comparison_source.tsv", sep="\t", index=False)

    _render_scatter_grid(dnase, statistics, "DNASE", "A", run_dir / "figures/FigureS6A.svg")
    _render_scatter_grid(atac, statistics, "ATAC", "B", run_dir / "figures/FigureS6B.svg")
    _render_paired_summary(statistics, run_dir / "figures/FigureS6C.svg")


def _element_statistics(merged: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (element, modality, model, model_version), group in merged.groupby(
        ["element", "modality", "model", "model_version"], sort=False, dropna=False
    ):
        clean = group[["mpra_effect", "ag_accessibility_effect"]].dropna()
        rho = stats.spearmanr(clean.mpra_effect, clean.ag_accessibility_effect)
        pearson = stats.pearsonr(clean.mpra_effect, clean.ag_accessibility_effect)
        z = np.arctanh(np.clip(float(rho.statistic), -.999999, .999999))
        se = 1 / np.sqrt(len(clean) - 3)
        rows.append({
            "element": element, "n": len(clean), "spearman_rho": float(rho.statistic),
            "spearman_p": float(rho.pvalue), "spearman_ci_low": float(np.tanh(z - 1.96 * se)),
            "spearman_ci_high": float(np.tanh(z + 1.96 * se)), "pearson_r": float(pearson.statistic),
            "pearson_p": float(pearson.pvalue),
            "direction_agreement": float((np.sign(clean.mpra_effect) == np.sign(clean.ag_accessibility_effect)).mean()),
            "modality": modality, "cell_line": ELEMENTS[element][3], "model": model, "model_version": model_version,
        })
    return pd.DataFrame(rows)


def _render_scatter_grid(scores: pd.DataFrame, statistics: pd.DataFrame, modality: str, letter: str, path: Path) -> None:
    available = [e for e in OTHER_ELEMENTS if not scores[scores.element.eq(e)].empty]
    fig, axes = plt.subplots(2, len(available), figsize=(2.2 * len(available), 4.2), squeeze=False)
    for row, model in enumerate(("ALL_FOLDS", "matched held-out")):
        for col, element in enumerate(available):
            ax = axes[row, col]
            panel = scores[scores.element.eq(element) & scores.model.eq(model)]
            ax.scatter(panel.mpra_effect, panel.ag_accessibility_effect, s=4, alpha=0.2, color=MODEL_COLORS[model], edgecolors="none")
            x, y = panel.mpra_effect.to_numpy(float), panel.ag_accessibility_effect.to_numpy(float)
            if len(x) > 1:
                slope, intercept = np.polyfit(x, y, 1)
                grid = np.linspace(x.min(), x.max(), 100)
                ax.plot(grid, slope * grid + intercept, color=MODEL_COLORS[model], lw=1.0)
            stat = statistics[statistics.element.eq(element) & statistics.model.eq(model) & statistics.modality.eq(modality)].iloc[0]
            label = f"ρ={stat.spearman_rho:.2f}"
            if model == "matched held-out":
                label += f"\n{HELDOUT_FOLD[element]}"
            ax.text(0.04, 0.96, label, transform=ax.transAxes, ha="left", va="top", fontsize=6.5)
            if row == 0:
                ax.set_title(DISPLAY[element], fontsize=8)
            if row == 1:
                ax.set_xlabel("Kircher MPRA effect", fontsize=6.5)
            if col == 0:
                ax.set_ylabel(f"{model}\nAG accessibility", fontsize=6.5)
    fig.suptitle(f"Figure S6{letter}: {modality}", fontsize=9)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)


def _render_paired_summary(statistics: pd.DataFrame, path: Path) -> None:
    ordered = []
    for element in OTHER_ELEMENTS:
        for modality in ("DNASE", "ATAC"):
            if not statistics[statistics.element.eq(element) & statistics.modality.eq(modality)].empty:
                ordered.append((element, modality))
    y = np.arange(len(ordered))[::-1]
    fig, ax = plt.subplots(figsize=(4.2, 3.0))
    for yi, (element, modality) in zip(y, ordered, strict=True):
        rows = statistics[statistics.element.eq(element) & statistics.modality.eq(modality)].set_index("model")
        ax.plot(
            [rows.loc["ALL_FOLDS", "spearman_rho"], rows.loc["matched held-out", "spearman_rho"]],
            [yi, yi], color="#B8B8B8", lw=0.9, zorder=1,
        )
        for model in ("ALL_FOLDS", "matched held-out"):
            row = rows.loc[model]
            ax.hlines(yi, row.spearman_ci_low, row.spearman_ci_high, color=MODEL_COLORS[model], lw=1.0, zorder=2)
            ax.scatter(row.spearman_rho, yi, s=22, color=MODEL_COLORS[model], edgecolor="white", lw=0.4, zorder=3)
    ax.axvline(0, color="#666666", ls="--", lw=0.65)
    ax.set_yticks(y, [f"{DISPLAY[e]} {'DNase' if m == 'DNASE' else 'ATAC'}" for e, m in ordered])
    ax.set_xlabel("Spearman rho: MPRA vs predicted accessibility")
    handles = [Line2D([0], [0], marker="o", ls="none", color=MODEL_COLORS[m], markersize=4.5, label=m) for m in MODEL_COLORS]
    fig.legend(handles=handles, frameon=False, loc="upper center", ncol=2, fontsize=6.5)
    fig.suptitle("Figure S6C", fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)
