#!/usr/bin/env python3
"""Render source-backed Figure 2 panels from frozen release tables.

Panels A and D are editable author-layout schematics in ``figures/assembled``.
Panels B, C, E and F are deterministic plots and require no AlphaGenome API
calls or public-data downloads.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["svg.hashsalt"] = "sort1-reanalysis-v1"
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[1]
SOURCE_DATA = REPO_ROOT / "outputs" / "source_data"
DEFAULT_OUTPUT = REPO_ROOT / "figures" / "rendered"
GENES = ("SORT1", "PSRC1", "CELSR2")


def _style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 8,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })


def _save(fig: plt.Figure, stem: str, output_dir: Path, **kwargs) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    svg = output_dir / f"{stem}.svg"
    svg_kwargs = dict(kwargs)
    svg_kwargs["metadata"] = {"Date": "2026-08-01", "Creator": "SORT1-reanalysis"}
    fig.savefig(svg, **svg_kwargs)
    fixed = datetime(2026, 8, 1, tzinfo=timezone.utc)
    pdf_kwargs = dict(kwargs)
    pdf_kwargs["metadata"] = {
        "Creator": "SORT1-reanalysis", "CreationDate": fixed, "ModDate": fixed,
    }
    fig.savefig(output_dir / f"{stem}.pdf", **pdf_kwargs)
    png_kwargs = dict(kwargs)
    png_kwargs["dpi"] = 600
    fig.savefig(output_dir / f"{stem}.png", **png_kwargs)
    plt.close(fig)
    return svg


def _centered(values: np.ndarray) -> TwoSlopeNorm:
    finite = values[np.isfinite(values)]
    maximum = float(np.max(np.abs(finite))) if finite.size else 1.0
    return TwoSlopeNorm(vmin=-maximum if maximum else -1, vcenter=0, vmax=maximum if maximum else 1)


def render_fig2b(output_dir: Path) -> Path:
    source = SOURCE_DATA / "Figure2B_top50_repair_outcomes"
    matrix = pd.read_csv(source / "matrix.csv").set_index("gene_symbol").loc[list(GENES)]
    columns = pd.read_csv(source / "columns.csv")
    if matrix.columns.tolist() != columns.column_label.tolist():
        raise ValueError("Figure 2B matrix and column manifest are out of order")
    values = matrix.to_numpy(float)
    width = max(8.7, 0.18 * matrix.shape[1] + 1.9)
    fig, ax = plt.subplots(figsize=(width, 2.65))
    fig.subplots_adjust(left=0.09, right=0.91, top=0.78, bottom=0.43)
    image = ax.pcolormesh(
        np.arange(matrix.shape[1] + 1), np.arange(matrix.shape[0] + 1), values,
        cmap="RdBu_r", norm=_centered(values), edgecolors="white",
        linewidth=0.35, antialiased=False, shading="flat",
    )
    ax.set_aspect("equal")
    ax.set_xlim(0, matrix.shape[1])
    ax.set_ylim(matrix.shape[0], 0)
    ax.set_title("Top repair outcomes", fontsize=10, pad=5)
    ax.set_yticks(np.arange(matrix.shape[0]) + 0.5, matrix.index, fontsize=8.5)
    ax.set_xticks(np.arange(matrix.shape[1]) + 0.5, matrix.columns,
                  rotation=90, ha="center", va="top", fontsize=6.2)
    ax.tick_params(axis="both", length=0, pad=2)
    for spine in ax.spines.values():
        spine.set_visible(False)
    colorbar = fig.colorbar(image, ax=ax, shrink=0.84, pad=0.012)
    colorbar.set_label("RNA change vs minor (%)", fontsize=8.5)
    colorbar.ax.tick_params(labelsize=7, length=2)
    return _save(fig, "Figure2B", output_dir, bbox_inches="tight")


def render_fig2c(output_dir: Path) -> Path:
    data = pd.read_csv(SOURCE_DATA / "Figure2C_deletion_grid.csv")
    data = data[
        data.geometry.eq("full_xy_grid")
        & data.target_index.eq(4)
        & data.gene_symbol.isin(GENES)
    ].copy()
    data["display"] = 100 * data["percent_change_vs_minor"]
    normalizer = _centered(data.display.to_numpy(float))
    fig, axes = plt.subplots(1, 3, figsize=(10.6, 3.7), constrained_layout=True)
    image = None
    for ax, gene in zip(axes, GENES, strict=True):
        matrix = data[data.gene_symbol.eq(gene)].pivot(
            index="upstream_bases", columns="downstream_bases", values="display"
        ).reindex(index=range(11), columns=range(11))
        if matrix.isna().any().any():
            raise ValueError(f"Incomplete Figure 2C deletion grid for {gene}")
        image = ax.imshow(matrix.to_numpy(), aspect="equal", origin="lower",
                          cmap="RdBu_r", norm=normalizer, extent=(-0.5, 10.5, -0.5, 10.5))
        ax.set_title(gene, fontsize=11)
        ax.set_xticks(range(11))
        ax.set_yticks(range(11))
        ax.set_xlabel("bases deleted downstream of cut site", fontsize=9.5)
        ax.tick_params(labelsize=8)
    axes[0].set_ylabel("bases deleted upstream of cut site", fontsize=9.5)
    colorbar = fig.colorbar(image, ax=axes, shrink=0.86, pad=0.02)
    colorbar.set_label("RNA change vs minor (%)", fontsize=9.5)
    colorbar.ax.tick_params(labelsize=8)
    return _save(fig, "Figure2C", output_dir, bbox_inches="tight")


def render_fig2e(output_dir: Path) -> Path:
    data = pd.read_csv(SOURCE_DATA / "Figure2E_kircher_correlations.tsv", sep="\t")
    data["overlaps_minor_motif"] = data.Position.between(109_274_967, 109_274_976)
    panels = [
        ("RNA", "ag_rna_3gene_mean_percent_change", r"AG $\Delta$ RNA %", (-31, 26), [-20, 0, 20]),
        ("ATAC", "ag_atac_mean_score", r"AG $\Delta$ ATAC", (-1.2, 0.45), [-1, 0, 0.4]),
        ("H3K27ac", "ag_h3k27ac_mean_score", r"AG $\Delta$ H3K27ac", (-2.1, 0.55), [-2, -1, 0]),
    ]
    x_column = "kircher_primary_log2_effect"
    fig = plt.figure(figsize=(3.0, 2.1), dpi=100)
    axes = [fig.add_axes([x0, 0.365, 0.72 / 3.0, 0.72 / 2.1]) for x0 in [0.105, 0.425, 0.745]]
    for ax, (title, y_column, ylabel, ylim, yticks) in zip(axes, panels, strict=True):
        subset = data.dropna(subset=[x_column, y_column])
        background = subset[~subset.overlaps_minor_motif]
        motif = subset[subset.overlaps_minor_motif]
        ax.scatter(background[x_column], background[y_column], s=3.2,
                   color="#809EEB", edgecolors="none", alpha=0.28, zorder=1)
        ax.scatter(motif[x_column], motif[y_column], s=4.6,
                   color="#C93636", edgecolors="none", alpha=0.82, zorder=3)
        x = subset[x_column].to_numpy(float)
        y = subset[y_column].to_numpy(float)
        slope, intercept = np.polyfit(x, y, 1)
        line_x = np.linspace(-3.2, 2.25, 100)
        ax.plot(line_x, slope * line_x + intercept, color="#222222", linewidth=0.85, zorder=2)
        correlation = stats.pearsonr(x, y).statistic
        ax.axhline(0, color="#666666", linewidth=0.45, zorder=0)
        ax.axvline(0, color="#666666", linewidth=0.45, zorder=0)
        ax.set_xlim(-3.15, 2.25)
        ax.set_ylim(*ylim)
        ax.set_xticks([-2, 0, 2])
        ax.set_yticks(yticks)
        ax.set_ylabel(ylabel, fontsize=4.8, labelpad=3.0)
        ax.set_title(f"{title}  r={correlation:.2f}", fontsize=5.8, pad=1.5)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(axis="both", labelsize=5.0, width=0.5, length=2.0, pad=1.2)
    fig.text(0.5, 0.215, "Kircher log2 effect", ha="center", va="center", fontsize=6.2)
    return _save(fig, "Figure2E", output_dir)


def render_fig2f(output_dir: Path) -> Path:
    source = SOURCE_DATA / "Figure2F_kircher_multielement"
    table = pd.read_csv(source / "kircher_multielement_element_statistics.tsv", sep="\t")
    summary = json.loads((source / "RESULTS.json").read_text())
    modality_colors = {"ATAC": "#1b9e77", "DNASE": "#7570b3"}
    markers = {"HepG2": "o", "K562": "s", "HeLa-S3": "^", "LNCaP clone FGC": "v"}
    labels = {"F9": "F9", "FOXE1": "FOXE1", "LDLR": "LDLR",
              "MYC_rs11986220": "MYC", "PKLR": "PKLR", "SORT1": "SORT1"}
    order = table[table.modality.eq("DNASE")].sort_values("spearman_rho").element.tolist()
    ymap = {element: index + 2 for index, element in enumerate(order)}
    fig, ax = plt.subplots(figsize=(5.8, 3.25))
    fig.subplots_adjust(left=0.27, right=0.75, bottom=0.19, top=0.98)
    for row in table.itertuples(index=False):
        y = ymap[row.element] + {"DNASE": -0.11, "ATAC": 0.11}[row.modality]
        color = modality_colors[row.modality]
        ax.hlines(y, row.spearman_ci_low, row.spearman_ci_high,
                  color=color, lw=1.25, alpha=0.92, zorder=1)
        ax.scatter(row.spearman_rho, y, color=color, marker=markers[row.cell_line],
                   s=32, edgecolor="white", linewidth=0.45, zorder=3)
    pooled_y = 0.65
    for modality, offset in (("DNASE", -0.09), ("ATAC", 0.09)):
        pooled = summary["matched_four_random_effects_by_modality"][modality]
        color = modality_colors[modality]
        ax.hlines(pooled_y + offset, pooled["pooled_95ci_low"], pooled["pooled_95ci_high"],
                  color=color, lw=1.6, zorder=1)
        ax.scatter(pooled["pooled_spearman_rho"], pooled_y + offset, color=color,
                   marker="D", s=38, edgecolor="white", linewidth=0.45, zorder=3)
    ax.axhline(1.45, color="0.82", lw=0.65)
    ax.axvline(0, color="0.35", lw=0.7, ls="--")
    ax.set_yticks([pooled_y] + [ymap[e] for e in order],
                  ["Matched subset (k=4)"] + [labels[e] for e in order])
    ax.set_xlabel("Spearman ρ: MPRA vs predicted accessibility")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="both", labelsize=7)
    ax.set_ylim(-0.35, max(ymap.values()) + 0.55)
    modality_handles = [Line2D([0], [0], color=modality_colors[m], lw=2,
                               label="DNase" if m == "DNASE" else "ATAC")
                        for m in ("DNASE", "ATAC")]
    first = ax.legend(handles=modality_handles, title="Output", loc="upper left",
                      bbox_to_anchor=(1.02, 0.98), frameon=False, fontsize=6.6,
                      title_fontsize=6.8, handlelength=1.7)
    ax.add_artist(first)
    cell_handles = [Line2D([0], [0], color="0.35", marker=marker, linestyle="none",
                           markersize=5, label=cell.replace(" clone FGC", ""))
                    for cell, marker in markers.items()]
    ax.legend(handles=cell_handles, title="Cell line", loc="upper left",
              bbox_to_anchor=(1.02, 0.57), frameon=False, fontsize=6.6,
              title_fontsize=6.8, handletextpad=0.6)
    return _save(fig, "Figure2F", output_dir, bbox_inches="tight")


RENDERERS = {"B": render_fig2b, "C": render_fig2c, "E": render_fig2e, "F": render_fig2f}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", choices=sorted(RENDERERS))
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    _style()
    for panel in ([args.panel] if args.panel else sorted(RENDERERS)):
        print(f"Figure 2{panel}: {RENDERERS[panel](args.output_dir)}")


if __name__ == "__main__":
    main()
