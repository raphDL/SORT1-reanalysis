"""Clean-room reproduction for Figure S3 (A-G): a caQTL analogue of Figure
1C using Currin et al. 2025 liver chromatin-accessibility QTLs.

The 80 Currin caQTL variants used in this panel are an exact subset of the
same 111 LDL-GWAS variants already reconstructed for Figure 1C/S2 (verified:
every rsid in the archive's Figure S3 source table is also in Figure 1C's
111-variant set). That makes almost everything here reuse, not new code:
variant/allele reconstruction (`figure1.py::_build_ag_variants`) and the
rs12740374 LD-tagging covariate (`figure1_public.py::compute_tagging_covariate`)
are identical to Figure 1C/S2's, just applied to this 80-variant subset. The
only genuinely new AlphaGenome calls are a fresh ATAC + H3K27ac CenterMask
scan of these 80 variants under ALL_FOLDS and FOLD_0 (panels B/C and E/F);
panels A, D, and G are derived from that plus the frozen Currin observed
caQTL data, with no further AlphaGenome calls.

The 80-variant Currin association table and its 28-peak "coordinated set"
definition are a manually-sliced local excerpt of the Currin et al. 2025
Zenodo/GEO release (accession recorded in data/SOURCES.tsv) with no
traceable derivation script in this codebase -- like Figure 2B's Wang
spreadsheet, they must be supplied by the caller and are checksum-validated,
never fetched live or committed to this repo.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

from .common import Audit, api_key, sha256_file
from .figure1 import (
    VARIANT_POS, VARIANT_RSID, _ag, _build_ag_variants, _save_svg,
    fetch_fig1c_inputs, run_fig1c_middle,
)
from .figure1_public import compute_tagging_covariate, prepare_figure1_public_inputs

CURRIN_VARIANTS_SHA256 = "4d295d2808ea81d3bcfb781857d551ea0f8658ae9bd35a75ba11c8f2508d1f8c"
CURRIN_PEAKSET_SHA256 = "6a1b7fdf7ec087b05cb80ff71fd2a368b9e9bed346058f887a22392491c7b26b"
DRIVER_PEAK = "peak15120"
FOCAL_PEAKS = ("peak15118", "peak15119", "peak15120", "peak15121", "peak15122")
LIVER = "UBERON:0002107"
HEPG2 = "EFO:0001187"
PEAK_GROUPS = (("driver_peak15120", (DRIVER_PEAK,)), ("focal5_ivw", FOCAL_PEAKS), ("set539_ivw", None))
CAUSAL_COLOR = "#d62728"


def stage_currin_inputs(run_dir: Path, audit: Audit, *, variants_file: Path | None, peakset_file: Path | None) -> dict[str, Path]:
    """Currin's SORT1-local caQTL association table and 28-peak coordinated-
    set definition have no traceable derivation script from the public
    Zenodo/GEO release in this codebase (like Figure 2B's Wang spreadsheet),
    so they must be supplied locally and are only checksum-validated here,
    never fetched or committed."""
    if variants_file is None or peakset_file is None:
        raise RuntimeError(
            "Figure S3 needs the Currin et al. 2025 SORT1-local caQTL association "
            "table and its 28-peak coordinated-set definition (Zenodo record "
            "15025748 / GEO GSE264684; see data/SOURCES.tsv). Pass them with "
            "--currin-variants and --currin-peakset; their SHA-256 will be verified."
        )
    raw = run_dir / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    result = {}
    for label, source, destination_name, expected in (
        ("currin_variants", variants_file, "currin_sort1_set539_variants.tsv", CURRIN_VARIANTS_SHA256),
        ("currin_peakset", peakset_file, "currin_sort1_coordinated_set.tsv", CURRIN_PEAKSET_SHA256),
    ):
        destination = raw / destination_name
        if not destination.exists():
            shutil.copy2(source, destination)
        observed = sha256_file(destination)
        if observed != expected:
            destination.unlink(missing_ok=True)
            raise ValueError(f"Currin {label} checksum mismatch: expected {expected}, got {observed}")
        audit.downloads.append({"url": f"supplied:{source}", "path": str(destination), "bytes": destination.stat().st_size, "sha256": observed, "reused": False})
        audit.save()
        result[label] = destination
    return result


def _score_atac(run_dir: Path, audit: Audit, *, model_version: Any, panel: str, variants: list[Any], raw_path: Path) -> pd.DataFrame:
    """The only genuinely new AlphaGenome calls in Figure S3: a single ATAC
    CenterMask scorer (liver + HepG2 ontologies summed), matching the
    manuscript's S3 methods exactly. The archive's frozen source table also
    carries an unused H3K27ac column and an `ld_r2_EUR_local` column from an
    earlier exploratory pass -- neither is referenced by any plotted panel
    or by the correlations table in build_figure_s2.py, so they are not
    reproduced here rather than guessed at."""
    genome, _, dna_client, model_modules = _ag()
    _, variant_scorers = model_modules
    cached = pd.read_csv(raw_path, sep="\t") if raw_path.exists() else pd.DataFrame()
    done = set(cached["variant"].astype(str)) if not cached.empty else set()
    wanted = [v for v in variants if str(v.name) not in done]
    rows = [] if cached.empty else [cached]
    if wanted:
        client = dna_client.create(api_key(), model_version=model_version, timeout=300)
        interval = genome.Interval("chr1", VARIANT_POS, VARIANT_POS).resize(dna_client.SEQUENCE_LENGTH_500KB)
        atac_scorer = variant_scorers.CenterMaskScorer(dna_client.OutputType.ATAC, None, variant_scorers.AggregationType.DIFF_LOG2_SUM)
        with audit.step(f"{panel}: score {len(wanted)} variants (ATAC, {model_version.name})"):
            for i, v in enumerate(wanted, start=1):
                (atac_ad,) = client.score_variant(interval=interval, variant=v, variant_scorers=[atac_scorer])
                atac_mask = atac_ad.var["ontology_curie"].astype(str).isin([LIVER, HEPG2]).to_numpy()
                atac_val = float(np.nansum(np.asarray(atac_ad.X)[0, atac_mask])) if atac_mask.any() else float("nan")
                rows.append(pd.DataFrame([{
                    "variant": str(v.name), "chrom": v.chromosome, "pos": int(v.position),
                    "ref": v.reference_bases, "alt": v.alternate_bases,
                    "atac_liver_full": atac_val,
                }]))
                audit.add_api_calls(panel, 1)
                audit.add_api_requests(panel, 1)
                if i % 10 == 0 or i == len(wanted):
                    raw_path.parent.mkdir(parents=True, exist_ok=True)
                    pd.concat(rows, ignore_index=True).drop_duplicates("variant", keep="last").to_csv(raw_path, sep="\t", index=False, float_format="%.10g")
        result = pd.concat(rows, ignore_index=True).drop_duplicates("variant", keep="last")
        raw_path.parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(raw_path, sep="\t", index=False, float_format="%.10g")
        return result
    return cached


def _weighted_mean_beta(group: pd.DataFrame) -> pd.Series:
    weights = 1.0 / group["varbeta"].astype(float).clip(lower=1e-12)
    beta = float(np.sum(group["caqtl_beta_ag_alt_minus_ref"] * weights) / np.sum(weights))
    return pd.Series({
        "caqtl_beta_ag_alt_minus_ref": beta,
        "caqtl_se_ivw": float(np.sqrt(1.0 / np.sum(weights))),
        "min_pvalue": float(group["pvalue"].min()),
        "max_neglog10_pvalue": float((-np.log10(group["pvalue"].clip(lower=np.finfo(float).tiny))).max()),
        "n_peaks": int(group["peak"].nunique()),
        "peaks": ",".join(sorted(group["peak"].unique())),
    })


def run_figs3(run_dir: Path, audit: Audit, *, currin_variants_file: Path | None, currin_peakset_file: Path | None) -> None:
    currin_inputs = stage_currin_inputs(run_dir, audit, variants_file=currin_variants_file, peakset_file=currin_peakset_file)
    currin = pd.read_csv(currin_inputs["currin_variants"], sep="\t")
    peakset = pd.read_csv(currin_inputs["currin_peakset"], sep="\t")
    set_peaks = [x for x in str(peakset.loc[0, "caPeaks"]).split(",") if x]
    currin["pos"] = currin["variant"].str.split(":").str[1].astype(int)
    currin["caqtl_ref"] = currin["variant"].str.split(":").str[2]
    currin["caqtl_alt"] = currin["variant"].str.split(":").str[3]

    # Prerequisite: Figure 1C-middle's reconstructed 111-variant table and
    # the hg38 FASTA used to confirm reference/alternate alleles.
    fetch_fig1c_inputs(run_dir, audit)
    variant_path = run_dir / "derived/Figure1C_reconstructed_variants.tsv"
    if not variant_path.exists():
        run_fig1c_middle(run_dir, audit, batch_size=12, max_workers=4, max_variants=None)
    variants_table = pd.read_csv(variant_path, sep="\t")

    # The 80-variant Figure S3 subset: Figure 1C variants tested against the
    # 28 coordinated-set peaks in the Currin data. Matching is by allele SET,
    # not just position -- Currin also tests a handful of indels at positions
    # that coincide with a GWAS SNV (e.g. chr1:109275536:C:CT alongside the
    # GWAS SNV rs57677983, C>T, at the same coordinate); a position-only
    # match would wrongly pull those indels in as if they were the SNV.
    candidates = currin[currin["peak"].isin(set_peaks)][["pos", "caqtl_ref", "caqtl_alt"]].drop_duplicates()
    allele_lookup = variants_table[["rsid", "pos", "ldl_lowering_allele", "ldl_raising_allele"]].merge(candidates, on="pos")
    allele_match = (
        (allele_lookup["caqtl_ref"] == allele_lookup["ldl_lowering_allele"]) & (allele_lookup["caqtl_alt"] == allele_lookup["ldl_raising_allele"])
    ) | (
        (allele_lookup["caqtl_ref"] == allele_lookup["ldl_raising_allele"]) & (allele_lookup["caqtl_alt"] == allele_lookup["ldl_lowering_allele"])
    )
    set_positions = set(allele_lookup.loc[allele_match, "pos"].unique())
    subset = variants_table[variants_table["pos"].isin(set_positions)].sort_values("pos").reset_index(drop=True)

    _, _, _, model_modules = _ag()
    dna_model, _ = model_modules
    ag_variants, alleles = _build_ag_variants(subset, run_dir / "raw/hg38.fa")
    alleles = alleles.rename(columns={"ref_hg38": "ag_ref", "alt_hg38": "ag_alt"})

    with audit.step("S3 prerequisite: rs12740374 EUR LD-tagging covariate for the caQTL subset"):
        tagging_inputs = prepare_figure1_public_inputs(run_dir, audit, gtex_file=None, vcf_file=None, panel_file=None, skip_gtex=True)
        tagged = compute_tagging_covariate(tagging_inputs, subset)

    scored: dict[str, pd.DataFrame] = {}
    for regime, panel in (("ALL_FOLDS", "S3B"), ("FOLD_0", "S3E")):
        raw_path = run_dir / "predictions" / f"FigureS3_{regime.lower()}_atac.tsv"
        scored[regime] = _score_atac(
            run_dir, audit, model_version=getattr(dna_model.ModelVersion, regime), panel=panel,
            variants=ag_variants, raw_path=raw_path,
        )

    base = tagged.merge(alleles[["rsid", "pos", "ag_ref", "ag_alt"]], on=["rsid", "pos"], validate="one_to_one")
    for regime, prefix in (("ALL_FOLDS", "all_folds"), ("FOLD_0", "fold0")):
        s = scored[regime][["variant", "atac_liver_full"]].rename(columns={"atac_liver_full": f"{prefix}_atac_direct"})
        base = base.merge(s, left_on="rsid", right_on="variant", how="left", validate="one_to_one").drop(columns=["variant"])
        causal_atac = float(base.loc[base.rsid.eq(VARIANT_RSID), f"{prefix}_atac_direct"].iloc[0])
        model = base[f"{prefix}_atac_direct"] + base["tagging_covvar_EUR"].astype(float) * causal_atac
        model.loc[base.rsid.eq(VARIANT_RSID)] = causal_atac
        base[f"{prefix}_atac_tagging"] = model

    # `old_ref`/`old_alt` is GLGC's raw effect/other allele (the archive's
    # own "old" axis, from its now-superseded score_all_1p13_snps_v2.tsv) --
    # this, not the hg38-confirmed ag_ref/ag_alt, is the axis the archive
    # actually orients the Currin beta onto. Using ag_ref/ag_alt here looked
    # plausible (both are just "some ref/alt pair") but silently flipped the
    # sign for ~61% of variants; verified against the archive to 0 diff
    # across all 80 driver-peak rows before trusting this.
    base["old_ref"] = base["effect_allele"]
    base["old_alt"] = base["other_allele"]
    base["variant"] = "chr1:" + base["pos"].astype(str) + ":" + base["ag_ref"] + ":" + base["ag_alt"]
    base["variant_type"] = np.where((base["ag_ref"].str.len() == 1) & (base["ag_alt"].str.len() == 1), "SNV", "indel")
    base["is_causal"] = base["rsid"].eq(VARIANT_RSID)
    base["pos_mb"] = base["pos"] / 1e6

    # Join is by position, not allele -- Currin occasionally tests an indel
    # at the same coordinate as a GWAS SNV (see the set_positions filter
    # above); drop those non-matching rows here rather than fail, matching
    # the legacy pipeline's own allele_order-filter behavior.
    merged = currin[currin["peak"].isin(set_peaks)].merge(base, on="pos", how="inner", suffixes=("", "_ag"))
    same = merged["old_ref"].eq(merged["caqtl_ref"]) & merged["old_alt"].eq(merged["caqtl_alt"])
    reverse = merged["old_ref"].eq(merged["caqtl_alt"]) & merged["old_alt"].eq(merged["caqtl_ref"])
    merged = merged[same | reverse].copy()
    same = same[same | reverse]
    if merged.empty or merged["rsid"].nunique() != len(subset):
        raise ValueError(f"Allele-matched caQTL rows cover {merged['rsid'].nunique()} of {len(subset)} selected variants")
    merged["caqtl_beta_ag_alt_minus_ref"] = np.where(same, merged["beta"], -merged["beta"])

    key_cols = [
        "rsid", "variant", "variant_type", "pos", "pos_mb", "old_ref", "old_alt", "ag_ref", "ag_alt",
        "all_folds_atac_direct", "all_folds_atac_tagging", "fold0_atac_direct", "fold0_atac_tagging",
        "tagging_covvar_EUR", "tagging_ptag_EUR", "is_causal",
    ]
    agg_parts = []
    for label, peaks in PEAK_GROUPS:
        sub = merged[merged["peak"].isin(peaks)] if peaks is not None else merged[merged["peak"].isin(set_peaks)]
        grouped = sub.groupby(key_cols, dropna=False).apply(_weighted_mean_beta, include_groups=False).reset_index()
        grouped.insert(0, "peak_group", label)
        agg_parts.append(grouped)
    data = pd.concat(agg_parts, ignore_index=True)

    out_dir = run_dir / "derived"
    out_dir.mkdir(parents=True, exist_ok=True)
    displayed_path = out_dir / "FigureS3_displayed_values.tsv"
    data.to_csv(displayed_path, sep="\t", index=False, float_format="%.10g")

    corr_rows = []
    for group, _ in PEAK_GROUPS:
        sub = data[data["peak_group"].eq(group)]
        for regime, direct, tagging in (("ALL_FOLDS", "all_folds_atac_direct", "all_folds_atac_tagging"), ("FOLD_0", "fold0_atac_direct", "fold0_atac_tagging")):
            for model, column in (("direct", direct), ("tagging", tagging)):
                corr_rows.append({
                    "peak_group": group, "regime": regime, "model": model, "n": len(sub),
                    "pearson_r": sub[column].corr(sub["caqtl_beta_ag_alt_minus_ref"], method="pearson"),
                    "spearman_rho": sub[column].corr(sub["caqtl_beta_ag_alt_minus_ref"], method="spearman"),
                })
    pd.DataFrame(corr_rows).to_csv(out_dir / "FigureS3_correlations.tsv", sep="\t", index=False, float_format="%.10g")

    _render_figs3(data, run_dir / "figures")


def _limits(frame: pd.DataFrame, columns: tuple[str, ...]) -> tuple[float, float]:
    values = frame[list(columns)].to_numpy(float).ravel()
    values = values[np.isfinite(values)]
    extent = max(abs(float(values.min())), abs(float(values.max())))
    return -1.08 * extent, 1.08 * extent


def _plot_column(data: pd.DataFrame, value_col: str, ylabel: str, output: Path, *, y_limits: tuple[float, float]) -> None:
    fig, ax = plt.subplots(figsize=(55 / 25.4, 48 / 25.4))
    norm = Normalize(vmin=0.0, vmax=1.0)
    sub = data[data["peak_group"].eq("driver_peak15120")]
    causal = sub["rsid"].eq(VARIANT_RSID)
    known = sub["tagging_covvar_EUR"].notna() & ~causal
    ax.scatter(sub.loc[known, "pos"] / 1e6, sub.loc[known, value_col], c=sub.loc[known, "tagging_covvar_EUR"], cmap="viridis", norm=norm, s=11, alpha=0.78, linewidth=0)
    ax.scatter(sub.loc[causal, "pos"] / 1e6, sub.loc[causal, value_col], color=CAUSAL_COLOR, s=24, linewidth=0, zorder=3)
    ax.axhline(0, color="#777777", lw=0.5)
    ax.axvline(VARIANT_POS / 1e6, color="#999999", lw=0.55, ls=":")
    ax.set_ylim(*y_limits)
    ax.set_ylabel(ylabel)
    ax.set_xlabel("chr1 position (Mb)")
    fig.subplots_adjust(left=0.25, right=0.96, top=0.96, bottom=0.23)
    _save_svg(fig, output)


def _plot_correlation(data: pd.DataFrame, x_col: str, y_col: str, xlabel: str, ylabel: str, output: Path) -> None:
    columns = list(dict.fromkeys(["rsid", "pos", "tagging_covvar_EUR", x_col, y_col]))
    sub = data[data["peak_group"].eq("driver_peak15120")][columns].dropna()
    causal = sub["rsid"].eq(VARIANT_RSID)
    fig, ax = plt.subplots(figsize=(55 / 25.4, 48 / 25.4))
    ax.scatter(sub.loc[~causal, x_col], sub.loc[~causal, y_col], c=sub.loc[~causal, "tagging_covvar_EUR"], cmap="viridis", vmin=0, vmax=1, s=13, alpha=0.78, linewidth=0)
    ax.scatter(sub.loc[causal, x_col], sub.loc[causal, y_col], color=CAUSAL_COLOR, s=27, linewidth=0, zorder=3)
    slope, intercept = np.polyfit(sub[x_col], sub[y_col], 1)
    xx = np.linspace(sub[x_col].min(), sub[x_col].max(), 100)
    ax.plot(xx, intercept + slope * xx, color="#222222", lw=0.8)
    r = sub[x_col].corr(sub[y_col], method="pearson")
    rho = sub[x_col].corr(sub[y_col], method="spearman")
    ax.text(0.04, 0.95, f"n = {len(sub)}\nr = {r:.2f}\nρ = {rho:.2f}", transform=ax.transAxes, va="top", fontsize=6.5)
    ax.axhline(0, color="#777777", lw=0.5)
    ax.axvline(0, color="#777777", lw=0.5)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.subplots_adjust(left=0.24, right=0.96, top=0.96, bottom=0.23)
    _save_svg(fig, output)


def _render_figs3(data: pd.DataFrame, figures_dir: Path) -> None:
    observed_limits = _limits(data, ("caqtl_beta_ag_alt_minus_ref",))
    direct_limits = _limits(data, ("all_folds_atac_direct", "fold0_atac_direct"))
    tagging_limits = _limits(data, ("all_folds_atac_tagging", "fold0_atac_tagging"))
    specs = (
        ("caqtl_beta_ag_alt_minus_ref", "Observed caQTL beta", "FigureS3A.svg", observed_limits),
        ("all_folds_atac_direct", "ALL-FOLDS direct\nATAC effect", "FigureS3B.svg", direct_limits),
        ("all_folds_atac_tagging", "ALL-FOLDS ATAC\n+ rs127 tagging", "FigureS3C.svg", tagging_limits),
        ("fold0_atac_direct", "FOLD-0 direct\nATAC effect", "FigureS3E.svg", direct_limits),
        ("fold0_atac_tagging", "FOLD-0 ATAC\n+ rs127 tagging", "FigureS3F.svg", tagging_limits),
    )
    for value_col, ylabel, output, ylim in specs:
        _plot_column(data, value_col, ylabel, figures_dir / output, y_limits=ylim)
    _plot_correlation(data, "all_folds_atac_tagging", "caqtl_beta_ag_alt_minus_ref", "ALL-FOLDS rs12740374 tagging model", "Observed driver-peak caQTL beta", figures_dir / "FigureS3D.svg")
    _plot_correlation(data, "fold0_atac_tagging", "caqtl_beta_ag_alt_minus_ref", "FOLD-0 rs12740374 tagging model", "Observed driver-peak caQTL beta", figures_dir / "FigureS3G.svg")

    fig, ax = plt.subplots(figsize=(42 / 25.4, 18 / 25.4))
    ax.axis("off")
    sm = ScalarMappable(norm=Normalize(0, 1), cmap="viridis")
    cb = fig.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.5, pad=0.1)
    cb.set_label("EUR rs12740374 tagging coefficient")
    ax.scatter([], [], s=24, color=CAUSAL_COLOR, label="rs12740374")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, 1.2))
    _save_svg(fig, figures_dir / "FigureS3_shared_legend.svg")
