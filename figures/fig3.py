#!/usr/bin/env python3
"""Render Figure 3 panels from compact source tables in outputs/source_data/.

Pure deterministic plotting: no AlphaGenome API calls, no public-data
downloads. Each panel function reproduces the exact plot code (figure size,
colors, line widths, axis limits, labels) from the working-archive script
named in MANIFEST.tsv, adapted only to read from this repository's compact
release tables instead of the legacy analysis-output directories.

Panels 3A, 3B, and 3C are NOT yet ported here: their legacy scripts recompute
or re-derive values from large or external inputs (a 171 MB / 18.5 MB
Zenodo-pending table for 3A, an 860-line combined analysis+plot script for
3B, and a live genome-FASTA/JASPAR PWM scan for 3C) rather than purely
replotting an already-frozen compact table. See MANIFEST.tsv and
MANIFEST_NOTES.md for their legacy script pointers.

Usage:
    python figures/fig3.py            # renders panels G, E, F
    python figures/fig3.py --panel G  # renders a single panel
"""

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
from matplotlib.lines import Line2D

REPO_ROOT = Path(__file__).resolve().parents[1]
SOURCE_DATA = REPO_ROOT / "outputs" / "source_data"
RENDERED = REPO_ROOT / "figures" / "rendered"


# ---------------------------------------------------------------------------
# Figure 3G: component-necessity audit (three-gene mean)
# Legacy: figure3_restructured/panel_F_component_necessity_expanded/make_panel.py
# ---------------------------------------------------------------------------

_G_BLUE = "#2478b5"
_G_ORDER = [
    "native", "full_315", "H1", "H2", "H3", "H4", "H5", "H6", "C/EBP",
    "upstream_arm", "downstream_arm", "both_arms",
]
_G_LABELS = {
    "native": "Intact", "full_315": "Full 315",
    "upstream_arm": "Up", "downstream_arm": "Down", "both_arms": "Both",
}


def render_fig3g() -> Path:
    table = pd.read_csv(SOURCE_DATA / "Figure3G_component_necessity.tsv", sep="\t")
    table = table.set_index("component")

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"],
        "font.size": 8.5, "axes.labelsize": 9, "xtick.labelsize": 7.2,
        "ytick.labelsize": 7.5, "axes.linewidth": 0.8,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })

    x = np.arange(len(_G_ORDER), dtype=float)
    fig, ax = plt.subplots(figsize=(87 / 25.4, 54 / 25.4))
    fig.subplots_adjust(left=0.20, right=0.99, bottom=0.43, top=0.74)
    means = np.array([table.loc[item, "mean_retention"] for item in _G_ORDER])
    sems = np.array([table.loc[item, "sem_retention"] for item in _G_ORDER])
    ax.errorbar(
        x, means, yerr=sems, fmt="o", color=_G_BLUE, markersize=4.4,
        linewidth=0, elinewidth=0.9, capsize=1.8, zorder=4,
    )

    ax.axhline(1, color="#777777", linestyle=":", linewidth=0.9)
    ax.axhline(0, color="#bbbbbb", linewidth=0.65)
    ax.set_xticks(x, [_G_LABELS.get(item, item) for item in _G_ORDER])
    ax.tick_params(axis="x", labelsize=6.5, pad=7)
    plt.setp(ax.get_xticklabels(), rotation=90, ha="center", va="top", rotation_mode="anchor")
    ax.set_ylabel("Mean RNA retention\n(3 genes)", fontsize=8.5, labelpad=3)
    ax.set_ylim(-0.04, 1.10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#dddddd", linewidth=0.5)
    ax.set_axisbelow(True)
    ax.text(0.5, 1.025, "Controls", transform=ax.get_xaxis_transform(), ha="center", fontsize=6.5)
    ax.text(5.0, 1.025, "Individual windows", transform=ax.get_xaxis_transform(), ha="center", fontsize=6.5)
    ax.text(10.0, 1.025, "Larger components", transform=ax.get_xaxis_transform(), ha="center", fontsize=6.5)

    return _save(fig, "Figure3G", bbox_inches="tight", pad_inches=0.10)


# ---------------------------------------------------------------------------
# Figure 3E: directional scramble recovery (three-gene mean + individual genes)
# Legacy: figure3_restructured/panel_E_directional_scramble_recovery/make_panel_E.py
# ---------------------------------------------------------------------------

_E_GENES = ["SORT1", "PSRC1", "CELSR2"]
_E_ARM_STYLE = {"Upstream": ("#2478b5", "-"), "Downstream": ("#d62728", "--")}


def render_fig3e() -> Path:
    table = pd.read_csv(SOURCE_DATA / "Figure3E_directional_recovery.tsv", sep="\t")

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica Neue", "Helvetica", "DejaVu Sans"],
        "font.size": 8, "axes.labelsize": 8.5, "xtick.labelsize": 7,
        "ytick.labelsize": 7, "axes.linewidth": 0.8, "svg.fonttype": "none",
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })

    fig, ax = plt.subplots(figsize=(65 / 25.4, 45 / 25.4))
    fig.subplots_adjust(left=0.19, right=0.98, bottom=0.22, top=0.97)

    for arm, (color, linestyle) in _E_ARM_STYLE.items():
        arm_table = table[table["arm"] == arm]
        for gene in _E_GENES:
            gene_curve = arm_table[arm_table["series"] == gene].sort_values("extent_bp")
            ax.plot(
                gene_curve["extent_bp"], gene_curve["mean_retention"],
                color=color, linestyle=linestyle, linewidth=0.7, alpha=0.35, zorder=2,
            )
        mean_curve = arm_table[arm_table["series"] == "3-gene mean"].sort_values("extent_bp")
        x = mean_curve["extent_bp"].to_numpy(float)
        y = mean_curve["mean_retention"].to_numpy(float)
        sem = mean_curve["sem_retention"].to_numpy(float)
        ax.plot(x, y, color=color, linestyle=linestyle, linewidth=2.2,
                label=f"{arm}, 3-gene mean", zorder=4)
        ax.fill_between(x, y - sem, y + sem, color=color, alpha=0.11, linewidth=0)

    ax.axhline(1.0, color="#777777", linewidth=0.85, linestyle=":")
    ax.axhline(0.0, color="#aaaaaa", linewidth=0.7)
    ax.set_xlim(0, 1000)
    ax.set_ylim(-0.04, 1.16)
    ax.set_xticks([0, 200, 400, 600, 800, 1000])
    ax.set_xticklabels(["0\nmotif only", "200", "400", "600", "800", "1000"])
    ax.set_xlabel("Native sequence restored (bp)")
    ax.set_ylabel("Mean liver RNA retention\nvs intact locus")
    ax.grid(axis="y", color="#dddddd", linewidth=0.5)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    handles = [
        Line2D([0], [0], color="#2478b5", linewidth=2.2, label="Upstream mean"),
        Line2D([0], [0], color="#d62728", linewidth=2.2, linestyle="--", label="Downstream mean"),
        Line2D([0], [0], color="#777777", linewidth=0.7, alpha=0.55, label="Individual genes"),
    ]
    ax.legend(handles=handles, frameon=False, fontsize=6.1, loc="upper right", handlelength=2.3)

    return _save(fig, "Figure3E", bbox_inches="tight", pad_inches=0.025)


# ---------------------------------------------------------------------------
# Figure 3F: SORT1-only frozen boundary scan (uniform 1-bp assessment)
# Legacy: figure3_restructured/panel_F_optimal_window_heatmap/make_uniform_wide_main_panel.py
# ---------------------------------------------------------------------------


def render_fig3f() -> Path:
    src = SOURCE_DATA / "Figure3F_boundary_grid"
    frame = pd.read_csv(src / "surface_summary_paired.csv")
    frame = frame[frame["gene"].eq("SORT1")].copy()
    upstream = np.arange(177, 191)
    downstream = np.arange(105, 161)
    matrix = (
        frame.pivot(index="upstream_bp", columns="downstream_bp", values="outside_mean_retention")
        .reindex(index=upstream, columns=downstream)
        .to_numpy()
    )
    if np.isnan(matrix).any():
        raise ValueError("The uniform one-base grid source data is incomplete.")

    selected_mean = pd.read_csv(src / "selected_mean_window.csv").iloc[0]
    selected_median = pd.read_csv(src / "selected_median_window.csv").iloc[0]

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica Neue", "DejaVu Sans"],
        "font.size": 8, "axes.labelsize": 9, "xtick.labelsize": 8,
        "ytick.labelsize": 8, "svg.fonttype": "none",
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })

    fig, ax = plt.subplots(figsize=(6.35, 2.05))
    fig.subplots_adjust(left=0.105, right=0.85, bottom=0.23, top=0.82)
    mesh = ax.pcolormesh(
        np.arange(104.5, 161.5), np.arange(176.5, 191.5), matrix,
        cmap="viridis", vmin=0.40, vmax=0.80, shading="flat", rasterized=False,
    )
    ax.contour(downstream, upstream, matrix, levels=[0.70], colors="white", linewidths=1.0)

    mean_length = int(selected_mean["window_length_bp"])
    median_length = int(selected_median["window_length_bp"])
    ax.scatter(
        [selected_mean["downstream_bp"]], [selected_mean["upstream_bp"]],
        marker="s", s=48, facecolors="none", edgecolors="#ef3b2c", linewidths=1.4,
        label=f"{mean_length} bp mean-only", zorder=5,
    )
    ax.scatter(
        [selected_median["downstream_bp"]], [selected_median["upstream_bp"]],
        marker="D", s=44, facecolors="none", edgecolors="#fdae6b", linewidths=1.4,
        label=f"{median_length} bp median-passing", zorder=5,
    )

    ax.set_xlim(104.5, 160.5)
    ax.set_ylim(176.5, 190.5)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Downstream extent retained (bp)")
    ax.set_ylabel("Upstream extent retained (bp)")
    ax.set_xticks([105, 115, 125, 135, 145, 155])
    ax.set_yticks([177, 180, 185, 190])
    ax.tick_params(width=0.7, length=3)
    for spine in ax.spines.values():
        spine.set_linewidth(0.7)

    ax.legend(
        loc="lower right", bbox_to_anchor=(0.995, 1.02), frameon=False, ncol=2,
        fontsize=7.5, handletextpad=0.35, columnspacing=1.0, borderaxespad=0,
    )
    ax.text(
        0.985, 0.07, "white line = 70%", transform=ax.transAxes,
        ha="right", va="bottom", color="white", fontsize=7.2,
    )

    cbar = fig.colorbar(mesh, ax=ax, fraction=0.035, pad=0.025)
    cbar.set_label("Mean SORT1 retention\nvs intact locus", labelpad=4)
    cbar.set_ticks([0.4, 0.5, 0.6, 0.7, 0.8])
    cbar.ax.tick_params(width=0.7, length=2.5)

    return _save(fig, "Figure3F", bbox_inches="tight", pad_inches=0.025)


def _save(fig, stem: str, **savefig_kwargs) -> Path:
    RENDERED.mkdir(parents=True, exist_ok=True)
    out = RENDERED / f"{stem}.svg"
    # Matplotlib otherwise embeds the wall-clock time in SVG/PDF metadata and
    # may generate process-dependent SVG object IDs.  Fix both so repeated
    # renders from identical source tables are byte-reproducible.
    svg_kwargs = dict(savefig_kwargs)
    svg_kwargs["metadata"] = {
        "Date": "2026-08-01",
        "Creator": "SORT1-reanalysis",
    }
    fig.savefig(out, **svg_kwargs)
    pdf_kwargs = dict(savefig_kwargs)
    fixed_time = datetime(2026, 8, 1, tzinfo=timezone.utc)
    pdf_kwargs["metadata"] = {
        "Creator": "SORT1-reanalysis",
        "CreationDate": fixed_time,
        "ModDate": fixed_time,
    }
    fig.savefig(RENDERED / f"{stem}.pdf", **pdf_kwargs)
    png_kwargs = dict(savefig_kwargs)
    png_kwargs["dpi"] = 600
    fig.savefig(RENDERED / f"{stem}.png", **png_kwargs)
    plt.close(fig)
    return out


RENDERERS = {"G": render_fig3g, "E": render_fig3e, "F": render_fig3f}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", choices=sorted(RENDERERS), default=None,
                         help="Render a single panel letter; default renders all ported panels.")
    args = parser.parse_args()

    panels = [args.panel] if args.panel else sorted(RENDERERS)
    for panel in panels:
        out = RENDERERS[panel]()
        print(f"Figure 3{panel}: {out}")


if __name__ == "__main__":
    main()
