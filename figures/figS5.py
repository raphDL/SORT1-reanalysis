#!/usr/bin/env python3
"""Render the split Figure S5 chromatin panels from frozen compact tables."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["svg.hashsalt"] = "sort1-reanalysis-v1"
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "outputs/source_data"
OUT = ROOT / "figures/rendered"
MODEL_ORDER = ("ALL_FOLDS", "FOLD_0")
MODEL_COLORS = {"ALL_FOLDS": "#377EB8", "FOLD_0": "#E68632"}


def render(panel: str) -> Path:
    if panel == "B":
        path, column, title = SOURCE / "FigureS5B_ATAC_substitutions.tsv", "ag_atac_mean_score", "ATAC"
    else:
        path, column, title = SOURCE / "FigureS5C_H3K27ac_substitutions.tsv", "ag_h3k27ac_mean_score", "H3K27ac"
    data = pd.read_csv(path, sep="\t")

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica Neue", "DejaVu Sans"],
        "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8.5,
        "xtick.labelsize": 7, "ytick.labelsize": 7,
        "axes.linewidth": 0.7, "svg.fonttype": "none",
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    fig, axes = plt.subplots(2, 1, figsize=(78 / 25.4, 88 / 25.4), sharex=True)
    fig.subplots_adjust(left=0.23, right=0.97, top=0.91, bottom=0.15, hspace=0.30)
    for ax, model in zip(axes, MODEL_ORDER):
        frame = data.loc[data["model"].eq(model)].copy()
        x = pd.to_numeric(frame["kircher_mean_log2_effect"], errors="coerce")
        y = pd.to_numeric(frame[column], errors="coerce")
        valid = x.notna() & y.notna()
        frame, x, y = frame.loc[valid], x.loc[valid], y.loc[valid]
        rs = frame["is_rs12740374_position"].astype(bool)
        motif = frame["in_minor_cebpa_motif"].astype(bool) & ~rs
        background = ~(motif | rs)
        ax.scatter(x[background], y[background], s=7, color="#8AA9E8", alpha=0.25,
                   edgecolors="none", rasterized=True)
        ax.scatter(x[motif], y[motif], s=11, color="#D94841", alpha=0.82,
                   edgecolors="none", rasterized=True, zorder=3)
        ax.scatter(x[rs], y[rs], s=28, marker="D", color="#111111",
                   edgecolors="white", linewidths=0.4, zorder=4)
        slope, intercept = np.polyfit(x, y, 1)
        xs = np.linspace(float(x.min()), float(x.max()), 150)
        ax.plot(xs, slope * xs + intercept, color=MODEL_COLORS[model], linewidth=1.2)
        pearson = stats.pearsonr(x, y).statistic
        ax.text(0.04, 0.95, f"r = {pearson:.2f}\nn = {len(x):,}", transform=ax.transAxes,
                ha="left", va="top", fontsize=7,
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8, "pad": 0.5})
        ax.text(-0.21, 0.5, model, transform=ax.transAxes, rotation=90,
                ha="center", va="center", color=MODEL_COLORS[model],
                fontweight="bold", fontsize=7.5)
        ax.axhline(0, color="#777777", linewidth=0.5)
        ax.axvline(0, color="#777777", linewidth=0.5)
        ax.grid(color="#E8E8E8", linewidth=0.45)
        ax.set_axisbelow(True)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_ylabel("AG HepG2 score")
    axes[0].set_title(title, fontweight="bold", pad=3)
    axes[-1].set_xlabel("Mean Kircher MPRA log2 effect")

    OUT.mkdir(parents=True, exist_ok=True)
    stem = OUT / f"FigureS5{panel}"
    svg_meta = {"Date": "2026-08-01", "Creator": "SORT1-reanalysis"}
    pdf_time = datetime(2026, 8, 1, tzinfo=timezone.utc)
    pdf_meta = {"Creator": "SORT1-reanalysis", "CreationDate": pdf_time, "ModDate": pdf_time}
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.025, metadata=svg_meta)
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.025, metadata=pdf_meta)
    fig.savefig(stem.with_suffix(".png"), bbox_inches="tight", pad_inches=0.025, dpi=600)
    plt.close(fig)
    return stem.with_suffix(".svg")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--panel", choices=("B", "C"))
    args = parser.parse_args()
    panels = (args.panel,) if args.panel else ("B", "C")
    for panel in panels:
        print(render(panel))


if __name__ == "__main__":
    main()
