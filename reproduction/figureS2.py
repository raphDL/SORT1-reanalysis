"""Clean-room reproduction for Figure S2 (A-C): held-out FOLD_0 replication
of the 111-variant liver eQTL / exon-mask RNA scan behind Figure 1C.

S2A is a re-plot of Figure 1C's already-computed ALL_FOLDS scan -- no new
AlphaGenome calls; it only reuses Figure 1C's own prerequisites (mirroring
how S1A/S1B reuse Figure 1D/1E). S2B is the only panel that makes fresh
AlphaGenome calls: the identical 111-variant, liver-restricted, exon-mask
RNA_SEQ scan as Figure 1C, but scored under the held-out FOLD_0 model
instead of ALL_FOLDS. S2C compares the two regimes' raw scores directly and
makes no new calls of its own.

The 111-variant scan itself (`score_gene_mask_variants`) and the GTEx/LD-
tagging merge (`attach_eqtl_and_tagging`) are shared, not duplicated, with
Figure 1C's own pipeline in figure1.py / figure1_public.py.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from .common import Audit
from .figure1 import GENES, VARIANT_RSID, _ag, _save_svg, run_fig1c_middle, score_gene_mask_variants
from .figure1_public import attach_eqtl_and_tagging, build_full_figure1c, prepare_figure1_public_inputs

CAUSAL_COLOR = "#d62728"
PLOT_GENE_ORDER = ("SORT1", "PSRC1", "CELSR2")  # matches the legacy panel layout, not GENES order


def _scatter_eqtl_panel(wide: pd.DataFrame, regime: str) -> tuple[plt.Figure, pd.DataFrame]:
    """The enlarged rs12740374 LD-tagging-model vs. GTEx liver eQTL scatter,
    shared by S2A (ALL_FOLDS) and S2B (FOLD_0)."""
    fig, axes = plt.subplots(1, 3, figsize=(150 / 25.4, 51 / 25.4))
    displayed = []
    for ax, gene in zip(axes, PLOT_GENE_ORDER, strict=True):
        # attach_eqtl_and_tagging only materializes the NaN-filled "for_plot"
        # variant of this column (Figure 1C's own plot needs the fallback).
        # After the dropna() below (which requires tagging_covvar_EUR to be
        # present), the fallback never fires, so the two columns are
        # numerically identical for every retained row; rename for parity
        # with the archive's column naming.
        src_col, xcol, ycol = f"ag_model_snp_plus_covvar_rs127_for_plot_{gene}", f"ag_model_snp_plus_covvar_rs127_{gene}", f"eqtl_liver_{gene}"
        sub = wide[["rsid", "pos", "tagging_covvar_EUR", src_col, ycol]].dropna().rename(columns={src_col: xcol}).copy()
        causal = sub["rsid"].eq(VARIANT_RSID)
        ax.scatter(
            sub.loc[~causal, xcol], sub.loc[~causal, ycol], c=sub.loc[~causal, "tagging_covvar_EUR"],
            cmap="viridis", vmin=0, vmax=1, alpha=0.82, s=14, linewidth=0,
        )
        ax.scatter(sub.loc[causal, xcol], sub.loc[causal, ycol], color=CAUSAL_COLOR, s=30, linewidth=0, zorder=3)
        r = pearsonr(sub[xcol], sub[ycol]).statistic
        rho = spearmanr(sub[xcol], sub[ycol]).statistic
        coef = np.polyfit(sub[xcol], sub[ycol], 1)
        xx = np.linspace(sub[xcol].min(), sub[xcol].max(), 100)
        ax.plot(xx, coef[0] * xx + coef[1], color="#222222", lw=0.9)
        ax.text(0.04, 0.94, f"n = {len(sub)}\nr = {r:.2f}\nρ = {rho:.2f}", transform=ax.transAxes, va="top", fontsize=6.5)
        ax.axhline(0, color="#777777", lw=0.5)
        ax.axvline(0, color="#777777", lw=0.5)
        ax.set_title(gene, fontweight="bold")
        ax.set_xlabel(f"{regime.replace('_', '-')} rs12740374\nLD-tagging model", fontsize=7)
        if ax is axes[0]:
            ax.set_ylabel("GTEx liver eQTL effect", fontsize=7)
        displayed.append(
            sub[["rsid", "pos", "tagging_covvar_EUR", xcol, ycol]].assign(
                gene=gene, model_regime=regime, pearson_r=r, spearman_rho=rho
            )
        )
    sm = plt.cm.ScalarMappable(norm=plt.Normalize(vmin=0, vmax=1), cmap="viridis")
    cb = fig.colorbar(sm, ax=axes, fraction=0.024, pad=0.025, aspect=22)
    cb.set_label("EUR rs12740374 tagging coefficient", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    fig.subplots_adjust(left=0.075, right=0.88, bottom=0.24, top=0.86, wspace=0.34)
    return fig, pd.concat(displayed, ignore_index=True)


def run_figs2a(run_dir: Path, audit: Audit) -> None:
    """S2A: enlarged ALL_FOLDS LD-tagging-model vs. GTEx eQTL scatter --
    reuses Figure 1C's own prerequisites; makes no new AlphaGenome calls."""
    wide_path = run_dir / "derived/Figure1C_eqtl_direct_tagging.tsv"
    if not wide_path.exists():
        middle_path = run_dir / "derived/Figure1C_middle_ag_scores.tsv"
        if not middle_path.exists():
            run_fig1c_middle(run_dir, audit, batch_size=12, max_workers=4, max_variants=None)
        with audit.step("S2A prerequisite: Figure 1C GTEx/1000G merge"):
            inputs = prepare_figure1_public_inputs(run_dir, audit, gtex_file=None, vcf_file=None, panel_file=None)
            build_full_figure1c(run_dir, inputs)
    with audit.step("S2A: enlarged ALL_FOLDS rs12740374 LD-tagging vs. GTEx liver eQTL"):
        wide = pd.read_csv(wide_path, sep="\t")
        fig, displayed = _scatter_eqtl_panel(wide, "ALL_FOLDS")
        out = run_dir / "derived/FigureS2A_gtex_vs_rs127_ld_tagging_all_folds.tsv"
        out.parent.mkdir(parents=True, exist_ok=True)
        displayed.to_csv(out, sep="\t", index=False)
        _save_svg(fig, run_dir / "figures/FigureS2A.svg")


def run_figs2b(run_dir: Path, audit: Audit) -> None:
    """S2B: the only S2 panel with new AlphaGenome calls -- reruns Figure
    1C's 111-variant scan under the held-out FOLD_0 model."""
    genome, _, dna_client, model_modules = _ag()
    dna_model, _ = model_modules
    fold0_middle_path = run_dir / "derived/FigureS2B_fold0_ag_scores.tsv"
    if not fold0_middle_path.exists():
        result, oriented = score_gene_mask_variants(
            run_dir, audit, model_version=dna_model.ModelVersion.FOLD_0, panel="S2B",
            raw_filename="FigureS2B_fold0_raw_track_scores.tsv",
            batch_size=12, max_workers=4, max_variants=None,
        )
        fold0_middle_path.parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(fold0_middle_path, sep="\t", index=False, float_format="%.10g")
        oriented.to_csv(run_dir / "derived/FigureS2B_fold0_oriented_track_scores.tsv", sep="\t", index=False, float_format="%.10g")
    wide_path = run_dir / "derived/FigureS2B_eqtl_direct_tagging.tsv"
    if not wide_path.exists():
        with audit.step("S2B: attach GTEx liver eQTL + LD-tagging covariate to FOLD_0 scores"):
            inputs = prepare_figure1_public_inputs(run_dir, audit, gtex_file=None, vcf_file=None, panel_file=None)
            fold0_base = pd.read_csv(fold0_middle_path, sep="\t")
            wide = attach_eqtl_and_tagging(run_dir, inputs, fold0_base)
            wide.to_csv(wide_path, sep="\t", index=False, float_format="%.10g")
    with audit.step("S2B: enlarged FOLD_0 rs12740374 LD-tagging vs. GTEx liver eQTL"):
        wide = pd.read_csv(wide_path, sep="\t")
        fig, displayed = _scatter_eqtl_panel(wide, "FOLD_0")
        out = run_dir / "derived/FigureS2B_gtex_vs_rs127_ld_tagging_fold_0.tsv"
        out.parent.mkdir(parents=True, exist_ok=True)
        displayed.to_csv(out, sep="\t", index=False)
        _save_svg(fig, run_dir / "figures/FigureS2B.svg")


def run_figs2c(run_dir: Path, audit: Audit) -> None:
    """S2C: direct ALL_FOLDS-vs-FOLD_0 score comparison -- reuses S2A/S2B's
    already-scored tables; makes no new AlphaGenome calls."""
    all_folds_path = run_dir / "derived/Figure1C_middle_ag_scores.tsv"
    eqtl_path = run_dir / "derived/Figure1C_eqtl_direct_tagging.tsv"
    fold0_path = run_dir / "derived/FigureS2B_fold0_ag_scores.tsv"
    if not all_folds_path.exists() or not eqtl_path.exists():
        run_figs2a(run_dir, audit)
    if not fold0_path.exists():
        run_figs2b(run_dir, audit)
    with audit.step("S2C: ALL_FOLDS vs. FOLD_0 direct score comparison"):
        all_folds_wide = pd.read_csv(all_folds_path, sep="\t")
        fold0_wide = pd.read_csv(fold0_path, sep="\t")
        eqtl_wide = pd.read_csv(eqtl_path, sep="\t")

        def _melt(wide: pd.DataFrame, value_name: str) -> pd.DataFrame:
            parts = [
                wide[["rsid", "pos", f"ag_rna_liver_{gene}"]].rename(columns={f"ag_rna_liver_{gene}": value_name}).assign(gene=gene)
                for gene in GENES
            ]
            return pd.concat(parts, ignore_index=True)

        paired = _melt(all_folds_wide, "ALL_FOLDS").merge(
            _melt(fold0_wide, "FOLD_0")[["rsid", "gene", "FOLD_0"]], on=["rsid", "gene"], validate="one_to_one"
        )
        eqtl_long = pd.concat(
            [
                eqtl_wide[["rsid", f"eqtl_liver_{gene}"]].rename(columns={f"eqtl_liver_{gene}": "GTEx_liver_eQTL"}).assign(gene=gene)
                for gene in GENES
            ],
            ignore_index=True,
        )
        paired = paired.merge(eqtl_long, on=["rsid", "gene"], how="left")

        fig, axes = plt.subplots(1, 3, figsize=(132 / 25.4, 48 / 25.4))
        rows = []
        for ax, gene in zip(axes, GENES, strict=True):
            sub = paired[paired["gene"].eq(gene)].copy()
            causal = sub["rsid"].eq(VARIANT_RSID)
            ax.scatter(sub.loc[~causal, "ALL_FOLDS"], sub.loc[~causal, "FOLD_0"], color="#aeb4ba", alpha=0.7, s=10, linewidth=0)
            ax.scatter(sub.loc[causal, "ALL_FOLDS"], sub.loc[causal, "FOLD_0"], color=CAUSAL_COLOR, s=25, linewidth=0, zorder=3)
            rho = spearmanr(sub["ALL_FOLDS"], sub["FOLD_0"]).statistic
            ax.text(0.04, 0.94, f"ρ = {rho:.2f}", transform=ax.transAxes, va="top")
            ax.axhline(0, color="#777777", lw=0.5)
            ax.axvline(0, color="#777777", lw=0.5)
            ax.set_title(gene, fontweight="bold")
            ax.set_xlabel("ALL_FOLDS exon-mask lnFC", fontsize=7)
            if ax is axes[0]:
                ax.set_ylabel("FOLD_0 exon-mask lnFC", fontsize=7)
            rs = sub.loc[causal].iloc[0]
            rows.append({
                "gene": gene, "n": len(sub),
                "spearman_rho_all_folds_vs_fold0": rho,
                "rs12740374_all_folds": rs["ALL_FOLDS"],
                "rs12740374_fold0": rs["FOLD_0"],
                "rs12740374_fold0_absolute_rank": int(
                    sub["FOLD_0"].abs().rank(ascending=False, method="min").loc[causal].iloc[0]
                ),
            })
        fig.tight_layout(w_pad=1.1)
        paired.to_csv(run_dir / "derived/FigureS2C_all_folds_vs_fold0.tsv", sep="\t", index=False, float_format="%.10g")
        pd.DataFrame(rows).to_csv(run_dir / "derived/FigureS2C_summary.tsv", sep="\t", index=False, float_format="%.10g")
        _save_svg(fig, run_dir / "figures/FigureS2C.svg")
