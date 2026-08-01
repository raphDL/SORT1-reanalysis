#!/usr/bin/env python3
"""Render source-backed Figure 1 panels from frozen release tables.

Panel A is an editable author-layout schematic in ``figures/assembled``.
Panel B can be rendered after its one-time numerical export with
``analysis/export_fig1b_source.py``; the frozen archive did not contain those
arrays. Panels C--F require no AlphaGenome calls or network access.

Usage:
    python figures/fig1.py
    python figures/fig1.py --panel D --output-dir /tmp/fig1-check
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
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.ticker import FormatStrFormatter, MaxNLocator

REPO_ROOT = Path(__file__).resolve().parents[1]
SOURCE_DATA = REPO_ROOT / "outputs" / "source_data"
DEFAULT_OUTPUT = REPO_ROOT / "figures" / "rendered"

VARIANT_POSITION = 109_274_968
CONTACT_GENES = ("CELSR2", "PSRC1", "SORT1")
TSS = {"CELSR2": 109_249_538, "PSRC1": 109_283_186, "SORT1": 109_397_918}
GENE_COLORS = {"CELSR2": "#59A14F", "PSRC1": "#E15759", "SORT1": "#4E79A7"}


def _style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.6,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
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


def _symmetric_limits(values: pd.Series) -> tuple[float, float]:
    maximum = float(np.nanmax(np.abs(values.to_numpy(float))))
    return (-1.08 * maximum, 1.08 * maximum)


def render_fig1b(output_dir: Path) -> Path:
    source = SOURCE_DATA / "Figure1B_locus_tracks" / "tracks.npz"
    if not source.exists():
        raise FileNotFoundError(
            "Figure 1B source arrays are absent. Run the authorized one-time "
            "analysis/export_fig1b_source.py export first."
        )
    tracks = np.load(source)
    specifications = [
        ("cebpa", "CEBPA ChIP", None),
        ("dnase", "DNase", None),
        ("atac", "ATAC", None),
        ("rna_plus", "RNA (+)", (0, 0.3)),
        ("rna_minus", "RNA (−)", (0, 1.8)),
    ]
    fig, axes = plt.subplots(5, 1, figsize=(4.1, 5.1), sharex=True)
    fig.subplots_adjust(left=0.20, right=0.97, top=0.98, bottom=0.12, hspace=0.18)
    for ax, (key, ylabel, ylim) in zip(axes, specifications, strict=True):
        start = int(tracks[f"{key}_start"][0])
        resolution = int(tracks[f"{key}_resolution"][0])
        ref = tracks[f"{key}_ref"]
        alt = tracks[f"{key}_alt"]
        positions = start + (np.arange(len(ref)) + 0.5) * resolution
        keep = (positions >= 109_250_000) & (positions <= 109_320_000)
        ax.plot(positions[keep], ref[keep], color="#1f77b4", lw=0.85)
        ax.plot(positions[keep], alt[keep], color="#d62728", lw=0.85)
        ax.axvspan(109_274_967.5, 109_274_968.5, color="goldenrod", alpha=0.35, lw=0)
        ax.set_ylabel(ylabel, fontsize=7)
        ax.tick_params(labelsize=6, length=2)
        if ylim:
            ax.set_ylim(*ylim)
        ax.spines[["top", "right"]].set_visible(False)
    axes[-1].set_xlim(109_250_000, 109_320_000)
    axes[-1].set_xlabel("chr1:109,250,000–109,320,000 (70 kb)", fontsize=7)
    axes[0].legend(
        handles=[
            Line2D([0], [0], color="#d62728", lw=1.2, label="ALT"),
            Line2D([0], [0], color="#1f77b4", lw=1.2, label="REF"),
        ],
        loc="upper right", frameon=False, fontsize=6, ncol=2,
    )
    return _save(fig, "Figure1B", output_dir, bbox_inches="tight")


def render_fig1c(output_dir: Path) -> Path:
    data = pd.read_csv(SOURCE_DATA / "Figure1C_eqtl_direct_tagging.tsv", sep="\t")
    genes = ("SORT1", "CELSR2", "PSRC1")
    norm = Normalize(0, 1)
    fig = plt.figure(figsize=(141 / 25.4, 130 / 25.4))
    grid = fig.add_gridspec(
        3, 3, hspace=0.30, wspace=0.54,
        left=0.125, right=0.865, top=0.925, bottom=0.115,
    )
    xlim = (float(data.pos_mb.min()) - 0.005, float(data.pos_mb.max()) + 0.005)
    causal = data.rsid.eq("rs12740374")
    tag = data["tagging_covvar_EUR"].notna()

    for row, gene in enumerate(genes):
        columns = (
            f"eqtl_liver_{gene}",
            f"ag_rna_liver_{gene}",
            f"ag_model_snp_plus_covvar_rs127_for_plot_{gene}",
        )
        ylabels = (f"{gene}\neQTL beta", "Exon-mask\nln fold change", "Exon-mask lnFC\n+ LD-tagging")
        titles = ("Observed liver eQTL", "AG single variant", "AG + LD-tagging")
        for col, (value_column, ylabel, title) in enumerate(zip(columns, ylabels, titles, strict=True)):
            ax = fig.add_subplot(grid[row, col])
            valid = data[value_column].notna()
            unknown = valid & ~tag & ~causal
            known = valid & tag & ~causal
            ax.scatter(data.loc[unknown, "pos_mb"], data.loc[unknown, value_column],
                       s=3.2, alpha=0.65, linewidths=0, color="#bdbdbd")
            ax.scatter(data.loc[known, "pos_mb"], data.loc[known, value_column],
                       s=3.2, alpha=0.82, linewidths=0,
                       c=data.loc[known, "tagging_covvar_EUR"], cmap="viridis", norm=norm)
            ax.scatter(data.loc[causal, "pos_mb"], data.loc[causal, value_column],
                       s=12, linewidths=0, color="#d62728", zorder=5)
            ax.axhline(0, color="black", lw=0.5)
            if col:
                ax.set_ylim(*_symmetric_limits(data.loc[valid, value_column]))
            else:
                values = data.loc[valid, value_column].to_numpy(float)
                pad = max(0.08 * float(np.ptp(values)), 1e-5)
                ax.set_ylim(min(0, float(values.min())) - pad, max(0, float(values.max())) + pad)
            ax.set_xlim(*xlim)
            ax.yaxis.set_major_locator(MaxNLocator(nbins=4, prune="both"))
            ax.xaxis.set_major_formatter(FormatStrFormatter("%.2f"))
            ax.xaxis.get_offset_text().set_visible(False)
            ax.tick_params(labelsize=6.0, pad=1.5)
            ax.set_ylabel(ylabel, fontsize=6.1, labelpad=1)
            if row == 0:
                ax.set_title(title, fontsize=7.1, pad=5)
            if row < 2:
                ax.set_xticklabels([])

    fig.text(0.495, 0.042, "chr1 position (Mb)", ha="center", va="center", fontsize=6.5)
    cax = fig.add_axes([0.885, 0.20, 0.018, 0.62])
    colorbar = fig.colorbar(ScalarMappable(norm=norm, cmap="viridis"), cax=cax)
    colorbar.set_label("EUR tagging\ncoef.", fontsize=5.8, labelpad=2)
    colorbar.ax.tick_params(labelsize=5.8, width=0.5, length=2)
    return _save(fig, "Figure1C", output_dir)


def _read_square_matrix(path: Path) -> tuple[np.ndarray, np.ndarray]:
    table = pd.read_csv(path, sep="\t")
    row_positions = table.iloc[:, 0].to_numpy(float)
    column_positions = np.asarray(table.columns[1:], dtype=float)
    if not np.allclose(row_positions, column_positions):
        raise ValueError(f"Row and column coordinates differ in {path.name}")
    values = table.iloc[:, 1:].to_numpy(float)
    if values.shape[0] != values.shape[1]:
        raise ValueError(f"Contact matrix is not square in {path.name}")
    return column_positions, values


def _contact_heatmap(path: Path, *, observed: bool, output_dir: Path, stem: str) -> Path:
    positions, values = _read_square_matrix(path)
    if observed:
        display = np.full(values.shape, np.nan, dtype=float)
        positive = values > 0
        display[positive] = np.log(values[positive])
    else:
        display = values
    fig = plt.figure(figsize=(60 / 25.4, 60 / 25.4))
    ax = fig.add_axes([0.16, 0.16, 0.66, 0.66])
    image = ax.imshow(
        display, origin="lower", aspect="equal", interpolation="nearest",
        cmap="inferno", norm=Normalize(-0.65, 2.2),
        extent=(positions[0], positions[-1], positions[0], positions[-1]),
    )
    guide = "cyan" if observed else "#222222"
    ax.axvline(0, color=guide, ls=":", lw=0.75)
    ax.axhline(0, color=guide, ls=":", lw=0.75)
    handles = []
    for gene in CONTACT_GENES:
        offset = (TSS[gene] - VARIANT_POSITION) / 1000
        for x0, y0 in ((offset, 0), (0, offset)):
            ax.add_patch(Rectangle((x0 - 2.6, y0 - 2.6), 5.2, 5.2,
                                   fill=False, edgecolor=GENE_COLORS[gene], linewidth=1.1))
        handles.append(Rectangle((0, 0), 1, 1, fill=False,
                                 edgecolor=GENE_COLORS[gene], linewidth=1.1, label=gene))
    ax.legend(handles=handles, frameon=False, fontsize=5.7, loc="lower center",
              bbox_to_anchor=(0.5, 1.015), ncol=3, columnspacing=0.8,
              handlelength=1.1, handletextpad=0.35)
    ax.set_xlabel("Position from rs12740374 (kb)", fontsize=6.8)
    ax.set_ylabel("Position from rs12740374 (kb)", fontsize=6.8)
    ax.tick_params(labelsize=6)
    cax = fig.add_axes([0.86, 0.25, 0.035, 0.48])
    cbar = fig.colorbar(image, cax=cax)
    cbar.set_ticks([-0.5, 0, 0.5, 1.0, 1.5, 2.0])
    cbar.ax.tick_params(labelsize=5.8)
    cbar.ax.set_title("log(O/E)", fontsize=5.6, pad=3)
    return _save(fig, stem, output_dir)


def render_fig1d(output_dir: Path) -> Path:
    return _contact_heatmap(SOURCE_DATA / "Figure1D_observed_hic.tsv",
                            observed=True, output_dir=output_dir, stem="Figure1D")


def render_fig1e(output_dir: Path) -> Path:
    return _contact_heatmap(SOURCE_DATA / "Figure1E_fold0_contact.tsv",
                            observed=False, output_dir=output_dir, stem="Figure1E")


def render_fig1f(output_dir: Path) -> Path:
    stats = pd.read_csv(SOURCE_DATA / "Figure1F_promoter_contact_percentiles.tsv", sep="\t")
    scales = ["Observed 4DN 2kb 3x3", "AlphaGenome FOLD_0 3x3"]
    pivot = stats[stats.scale.isin(scales)].pivot(
        index="gene", columns="scale", values="same_distance_percentile"
    ).loc[list(CONTACT_GENES)]
    fig = plt.figure(figsize=(60 / 25.4, 60 / 25.4))
    ax = fig.add_axes([0.25, 0.18, 0.62, 0.70])
    x = np.arange(3)
    width = 0.24
    bars = [
        ax.bar(x - width / 2, pivot[scales[0]], width, color="#555555", label="Experimental"),
        ax.bar(x + width / 2, pivot[scales[1]], width, color="#5B8CC0", label="FOLD_0 predicted"),
    ]
    ax.axhline(95, color="#B22222", ls="--", lw=0.7)
    ax.text(2.48, 96.3, "95th", color="#B22222", fontsize=5.4, ha="right", va="bottom")
    for series_index, collection in enumerate(bars):
        for gene_index, bar in enumerate(collection):
            value = float(bar.get_height())
            offset = 4.8 if series_index == 1 and gene_index == 2 else 2.0
            ax.text(bar.get_x() + bar.get_width() / 2, min(value + offset, 106),
                    f"{value:.1f}", ha="center", va="bottom", fontsize=5.5)
    ax.set_xticks(x, CONTACT_GENES)
    ax.set_ylim(0, 108)
    ax.set_ylabel("Matched-distance percentile", fontsize=6.8)
    ax.tick_params(labelsize=6)
    ax.legend(frameon=False, fontsize=5.8, loc="lower center",
              bbox_to_anchor=(0.5, 1.03), ncol=2, columnspacing=0.9, handlelength=1.3)
    return _save(fig, "Figure1F", output_dir)


RENDERERS = {"B": render_fig1b, "C": render_fig1c, "D": render_fig1d,
             "E": render_fig1e, "F": render_fig1f}
DEFAULT_PANELS = ("C", "D", "E", "F")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", choices=sorted(RENDERERS))
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    _style()
    for panel in ([args.panel] if args.panel else DEFAULT_PANELS):
        print(f"Figure 1{panel}: {RENDERERS[panel](args.output_dir)}")


if __name__ == "__main__":
    main()
