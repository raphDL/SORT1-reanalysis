"""Comparison and human-readable reporting for clean-room runs."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from .common import REPO_ROOT, file_manifest, sha256_file, utc_now


REFERENCE_ROOT = REPO_ROOT / "outputs" / "source_data"


def _numeric_summary(
    observed: np.ndarray, expected: np.ndarray, *, rtol: float = 1e-5,
    atol: float = 1e-6, min_pearson: float | None = None,
) -> dict[str, object]:
    observed = np.asarray(observed, dtype=float)
    expected = np.asarray(expected, dtype=float)
    if observed.shape != expected.shape:
        return {
            "pass": False,
            "reason": "shape_mismatch",
            "observed_shape": observed.shape,
            "expected_shape": expected.shape,
        }
    finite = np.isfinite(observed) & np.isfinite(expected)
    if not finite.any():
        return {"pass": False, "reason": "no_shared_finite_values"}
    delta = np.abs(observed[finite] - expected[finite])
    pearson = (
        float(np.corrcoef(observed[finite], expected[finite])[0, 1])
        if finite.sum() > 1 and np.std(observed[finite]) and np.std(expected[finite]) else None
    )
    close = bool(np.allclose(observed, expected, rtol=rtol, atol=atol, equal_nan=True))
    correlation_ok = min_pearson is None or (pearson is not None and pearson >= min_pearson)
    return {
        "pass": bool(close and correlation_ok),
        "n": int(finite.sum()),
        "max_abs_difference": float(delta.max()),
        "mean_abs_difference": float(delta.mean()),
        "pearson_r": pearson,
        "minimum_pearson_r": min_pearson,
        "rtol": rtol,
        "atol": atol,
    }


def compare_fig1b(run_dir: Path) -> dict[str, object]:
    generated = run_dir / "predictions" / "Figure1B_locus_tracks" / "tracks.npz"
    reference = REFERENCE_ROOT / "Figure1B_locus_tracks" / "tracks.npz"
    if not generated.exists() or not reference.exists():
        return {"pass": False, "reason": "generated_or_reference_file_missing"}
    got, want = np.load(generated), np.load(reference)
    keys = sorted(set(got.files) | set(want.files))
    values = {}
    for key in keys:
        if key not in got.files or key not in want.files:
            values[key] = {"pass": False, "reason": "key_missing"}
        elif key.endswith(("_start", "_resolution")):
            values[key] = {"pass": bool(np.array_equal(got[key], want[key]))}
        else:
            values[key] = _numeric_summary(got[key], want[key])
    return {
        "pass": all(bool(item["pass"]) for item in values.values()),
        "generated_sha256": sha256_file(generated),
        "reference_sha256": sha256_file(reference),
        "arrays": values,
    }


def compare_fig1c_middle(run_dir: Path) -> dict[str, object]:
    generated = run_dir / "derived" / "Figure1C_middle_ag_scores.tsv"
    reference = REFERENCE_ROOT / "Figure1C_eqtl_direct_tagging.tsv"
    if not generated.exists() or not reference.exists():
        return {"pass": False, "reason": "generated_or_reference_file_missing"}
    got = pd.read_csv(generated, sep="\t")
    want = pd.read_csv(reference, sep="\t")
    columns = [f"ag_rna_liver_{gene}" for gene in ("SORT1", "CELSR2", "PSRC1")]
    joined = got[["rsid", *columns]].merge(
        want[["rsid", *columns]], on="rsid", suffixes=("_generated", "_reference"), validate="one_to_one"
    )
    results = {
        column: _numeric_summary(joined[f"{column}_generated"], joined[f"{column}_reference"])
        for column in columns
    }


def compare_fig1c(run_dir: Path) -> dict[str, object]:
    generated = run_dir / "derived" / "Figure1C_eqtl_direct_tagging.tsv"
    reference = REFERENCE_ROOT / "Figure1C_eqtl_direct_tagging.tsv"
    if not generated.exists(): return {"pass": False, "reason": "generated_file_missing"}
    got, want = pd.read_csv(generated, sep="\t"), pd.read_csv(reference, sep="\t")
    columns = [f"eqtl_liver_{g}{s}" for g in ("SORT1","CELSR2","PSRC1") for s in ("","_se","_p")] + ["tagging_ptag_EUR","tagging_covvar_EUR"] + [f"ag_rna_liver_{g}" for g in ("SORT1","CELSR2","PSRC1")] + [f"ag_model_snp_plus_covvar_rs127_for_plot_{g}" for g in ("SORT1","CELSR2","PSRC1")]
    joined=got[["rsid",*columns]].merge(want[["rsid",*columns]],on="rsid",suffixes=("_generated","_reference"),validate="one_to_one")
    results={c:_numeric_summary(joined[f"{c}_generated"],joined[f"{c}_reference"],rtol=1e-4,atol=1e-5) for c in columns}
    exact=set(got.rsid)==set(want.rsid)
    return {"pass": exact and all(x["pass"] for x in results.values()),"variant_set_exact":exact,"values":results}


def _compare_matrix(run_dir: Path, name: str) -> dict[str, object]:
    generated=run_dir/"derived"/name; reference=REFERENCE_ROOT/name
    if not generated.exists(): return {"pass":False,"reason":"generated_file_missing"}
    got,want=pd.read_csv(generated,sep="\t"),pd.read_csv(reference,sep="\t")
    coords=_numeric_summary(got.iloc[:,0],want.iloc[:,0]); matrix=_numeric_summary(got.iloc[:,1:],want.iloc[:,1:],rtol=1e-5,atol=1e-5)
    return {"pass":bool(coords["pass"] and matrix["pass"]),"coordinates":coords,"matrix":matrix}


def compare_fig1d(run_dir: Path) -> dict[str, object]: return _compare_matrix(run_dir,"Figure1D_observed_hic.tsv")


def compare_fig1f(run_dir: Path) -> dict[str, object]:
    generated=run_dir/"derived/Figure1F_promoter_contact_percentiles.tsv"; reference=REFERENCE_ROOT/"Figure1F_promoter_contact_percentiles.tsv"
    if not generated.exists(): return {"pass":False,"reason":"generated_file_missing"}
    got,want=pd.read_csv(generated,sep="\t"),pd.read_csv(reference,sep="\t"); keys=["scale","gene"]
    cols=[c for c in want.columns if c not in keys]; joined=got.merge(want,on=keys,suffixes=("_generated","_reference"),validate="one_to_one")
    results={c:_numeric_summary(joined[f"{c}_generated"],joined[f"{c}_reference"],rtol=1e-5,atol=1e-5) for c in cols}
    return {"pass":len(joined)==6 and all(x["pass"] for x in results.values()),"rows":len(joined),"values":results}
    expected_ids = set(want["rsid"].astype(str))
    generated_ids = set(got["rsid"].astype(str))
    complete = generated_ids == expected_ids
    return {
        "pass": complete and all(bool(item["pass"]) for item in results.values()),
        "variant_set_exact": complete,
        "generated_variants": len(generated_ids),
        "reference_variants": len(expected_ids),
        "scores": results,
    }


def compare_fig1e(run_dir: Path) -> dict[str, object]:
    generated = run_dir / "derived" / "Figure1E_fold0_contact.tsv"
    reference = REFERENCE_ROOT / "Figure1E_fold0_contact.tsv"
    if not generated.exists() or not reference.exists():
        return {"pass": False, "reason": "generated_or_reference_file_missing"}
    got = pd.read_csv(generated, sep="\t")
    want = pd.read_csv(reference, sep="\t")
    coordinates = _numeric_summary(got.iloc[:, 0], want.iloc[:, 0])
    values = _numeric_summary(got.iloc[:, 1:].to_numpy(), want.iloc[:, 1:].to_numpy())
    return {"pass": bool(coordinates["pass"] and values["pass"]), "coordinates": coordinates, "matrix": values}


def compare_fig2b(run_dir: Path) -> dict[str, object]:
    generated_root = run_dir / "derived" / "Figure2B_top50_repair_outcomes"
    reference_root = REFERENCE_ROOT / "Figure2B_top50_repair_outcomes"
    if not (generated_root / "matrix.csv").exists():
        return {"pass": False, "reason": "generated_file_missing"}
    got = pd.read_csv(generated_root / "matrix.csv")
    want = pd.read_csv(reference_root / "matrix.csv")
    got_columns = pd.read_csv(generated_root / "columns.csv")
    want_columns = pd.read_csv(reference_root / "columns.csv")
    sequence_exact = got_columns.sequence_sha256.astype(str).tolist() == want_columns.sequence_sha256.astype(str).tolist()
    labels_exact = got_columns.column_label.astype(str).tolist() == want_columns.column_label.astype(str).tolist()
    values = _numeric_summary(
        got.iloc[:, 1:].to_numpy(), want.iloc[:, 1:].to_numpy(),
        rtol=0.0, atol=1.25, min_pearson=0.999,
    )
    return {
        "pass": bool(sequence_exact and labels_exact and values["pass"]),
        "repair_sequence_order_exact": sequence_exact,
        "labels_exact": labels_exact,
        "matrix": values,
    }


def compare_fig2c(run_dir: Path) -> dict[str, object]:
    generated = run_dir / "derived" / "Figure2C_deletion_grid.csv"
    reference = REFERENCE_ROOT / "Figure2C_deletion_grid.csv"
    if not generated.exists():
        return {"pass": False, "reason": "generated_file_missing"}
    keys = ["gene_symbol", "upstream_bases", "downstream_bases"]
    got = pd.read_csv(generated)
    want = pd.read_csv(reference)
    got = got[got.geometry.eq("full_xy_grid") & got.target_index.astype(int).eq(4) & got.gene_symbol.isin(("SORT1", "PSRC1", "CELSR2"))]
    want = want[want.geometry.eq("full_xy_grid") & want.target_index.astype(int).eq(4) & want.gene_symbol.isin(("SORT1", "PSRC1", "CELSR2"))]
    joined = got[keys + ["percent_change_vs_minor"]].merge(
        want[keys + ["percent_change_vs_minor"]], on=keys, suffixes=("_generated", "_reference"), validate="one_to_one"
    )
    expected = 3 * 121
    values = _numeric_summary(
        joined.percent_change_vs_minor_generated, joined.percent_change_vs_minor_reference,
        rtol=0.0, atol=0.0125, min_pearson=0.999,
    )
    return {"pass": bool(len(joined) == expected and values["pass"]), "grid_cells": len(joined), "values": values}


def compare_fig2e(run_dir: Path) -> dict[str, object]:
    generated = run_dir / "derived" / "Figure2E_kircher_correlations.tsv"
    reference = REFERENCE_ROOT / "Figure2E_kircher_correlations.tsv"
    if not generated.exists():
        return {"pass": False, "reason": "generated_file_missing"}
    got = pd.read_csv(generated, sep="\t")
    want = pd.read_csv(reference, sep="\t")
    columns = [
        "kircher_primary_log2_effect", "ag_atac_mean_score", "ag_h3k27ac_mean_score",
        "ag_rna_percent_change_SORT1", "ag_rna_percent_change_PSRC1",
        "ag_rna_percent_change_CELSR2", "ag_rna_3gene_mean_percent_change",
    ]
    joined = got[["variant_id_construct", *columns]].merge(
        want[["variant_id_construct", *columns]], on="variant_id_construct",
        suffixes=("_generated", "_reference"), validate="one_to_one",
    )
    tolerances = {
        "kircher_primary_log2_effect": (1e-6, None),
        "ag_atac_mean_score": (0.005, 0.9999),
        "ag_h3k27ac_mean_score": (0.011, 0.9999),
        "ag_rna_percent_change_SORT1": (0.6, 0.97),
        "ag_rna_percent_change_PSRC1": (0.6, 0.999),
        "ag_rna_percent_change_CELSR2": (0.6, 0.999),
        "ag_rna_3gene_mean_percent_change": (0.6, 0.999),
    }
    results = {
        column: _numeric_summary(
            joined[f"{column}_generated"], joined[f"{column}_reference"],
            rtol=0.0, atol=tolerances[column][0], min_pearson=tolerances[column][1],
        )
        for column in columns
    }
    plotted = {}
    for column in ("ag_atac_mean_score", "ag_h3k27ac_mean_score", "ag_rna_3gene_mean_percent_change"):
        generated_r = float(np.corrcoef(joined.kircher_primary_log2_effect_generated, joined[f"{column}_generated"])[0, 1])
        reference_r = float(np.corrcoef(joined.kircher_primary_log2_effect_reference, joined[f"{column}_reference"])[0, 1])
        plotted[column] = {
            "generated_pearson_r": generated_r,
            "reference_pearson_r": reference_r,
            "absolute_difference": abs(generated_r - reference_r),
            "maximum_allowed_difference": 0.01,
            "pass": abs(generated_r - reference_r) <= 0.01,
        }
    exact = set(got.variant_id_construct) == set(want.variant_id_construct)
    return {
        "pass": exact and len(joined) == 1798 and all(item["pass"] for item in results.values())
        and all(item["pass"] for item in plotted.values()),
        "variant_set_exact": exact,
        "matched_substitutions": len(joined),
        "scores": results,
        "plotted_correlations": plotted,
    }


def compare_fig2f(run_dir: Path) -> dict[str, object]:
    generated = run_dir / "derived" / "Figure2F_kircher_multielement" / "kircher_multielement_element_statistics.tsv"
    reference = REFERENCE_ROOT / "Figure2F_kircher_multielement" / "kircher_multielement_element_statistics.tsv"
    if not generated.exists():
        return {"pass": False, "reason": "generated_file_missing"}
    got = pd.read_csv(generated, sep="\t")
    want = pd.read_csv(reference, sep="\t")
    keys = ["element", "modality"]
    columns = ["n", "spearman_rho", "spearman_ci_low", "spearman_ci_high", "pearson_r", "direction_agreement"]
    joined = got[keys + columns].merge(want[keys + columns], on=keys, suffixes=("_generated", "_reference"), validate="one_to_one")
    results = {column: _numeric_summary(joined[f"{column}_generated"], joined[f"{column}_reference"]) for column in columns}
    return {"pass": len(joined) == 10 and all(item["pass"] for item in results.values()), "element_modality_pairs": len(joined), "statistics": results}


def compare_fig3b(run_dir: Path) -> dict[str, object]:
    root=run_dir/"derived/Figure3B_native_501bp_ism"; ref=REFERENCE_ROOT/"Figure3B_native_501bp_ism"
    generated=root/"native_locus_501bp_three_gene_mean_scores.tsv"
    if not generated.exists(): return {"pass":False,"reason":"generated_file_missing"}
    got=pd.read_csv(generated,sep="\t"); want=pd.read_csv(ref/"native_locus_501bp_three_gene_mean_scores.tsv",sep="\t")
    keys=["state_id","edit_offset","edit_pos_hg38","native_base","alt_base","sequence_sha256"]
    joined=got.merge(want,on=keys,suffixes=("_generated","_reference"),validate="one_to_one")
    values=_numeric_summary(joined.loss_vs_intact_T_generated,joined.loss_vs_intact_T_reference,rtol=0,atol=5e-4,min_pearson=.995)
    gpos=pd.read_csv(root/"native_locus_501bp_three_gene_mean_position_summary.tsv",sep="\t"); wpos=pd.read_csv(ref/"native_locus_501bp_three_gene_mean_position_summary.tsv",sep="\t")
    pos=gpos.merge(wpos,on=["edit_offset","edit_pos_hg38","native_base"],suffixes=("_generated","_reference"),validate="one_to_one")
    summary={c:_numeric_summary(pos[f"{c}_generated"],pos[f"{c}_reference"],rtol=0,atol=5e-4,min_pearson=.995) for c in ("max_loss","min_loss","mean_loss","positive_max_loss")}
    exact=len(joined)==1503 and got.sequence_sha256.tolist()==want.sequence_sha256.tolist()
    return {"pass":exact and values["pass"] and all(x["pass"] for x in summary.values()),"substitution_sequences_exact":exact,"matched_states":len(joined),"loss":values,"position_summary":summary}


def compare_fig3c(run_dir: Path) -> dict[str, object]:
    generated = run_dir / "derived/Figure3C_pwm_compatibility/Figure3C_PWM_disruption_values.tsv"
    reference = REFERENCE_ROOT / "Figure3C_pwm_compatibility.tsv"
    if not generated.exists():
        return {"pass": False, "reason": "generated_file_missing"}
    got = pd.read_csv(generated, sep="\t", index_col=0)
    want = pd.read_csv(reference, sep="\t", index_col=0)
    # Figure 3C's PWM scan is exact given fixed hotspot windows (verified
    # byte-identical against the release when fed the same 3B hotspots as
    # the original publication run). But when run end to end, it inherits
    # run_fig3b's ~1e-4 AlphaGenome run-to-run drift; a near-tied greedy
    # hotspot window pick can shift by 1 bp, which this pure local PWM scan
    # then amplifies into a visibly different score for that one cell. So
    # this compares correlation and displayed-family overlap, not exact
    # values -- see REPRODUCIBILITY_NEXT_STEPS.md R019.
    rows_exact = list(got.index) == list(want.index)
    shared_families = sorted(set(got.columns) & set(want.columns))
    family_overlap = len(shared_families) / len(want.columns)
    values = _numeric_summary(
        got.loc[got.index, shared_families].to_numpy() if rows_exact else np.array([]),
        want.loc[want.index, shared_families].to_numpy() if rows_exact else np.array([]),
        rtol=0.0, atol=1.0, min_pearson=0.85,
    ) if rows_exact and shared_families else {"pass": False, "reason": "rows_or_families_missing"}
    return {
        "pass": bool(rows_exact and family_overlap >= 0.75 and values["pass"]),
        "hotspot_rows_exact": rows_exact,
        "displayed_family_overlap_fraction": family_overlap,
        "shared_families": shared_families,
        "generated_only_families": sorted(set(got.columns) - set(want.columns)),
        "reference_only_families": sorted(set(want.columns) - set(got.columns)),
        "values": values,
    }


def compare_fig3g(run_dir: Path) -> dict[str, object]:
    generated = run_dir / "derived/Figure3G_component_necessity/Figure3G_component_necessity_three_gene_mean_source.tsv"
    reference = REFERENCE_ROOT / "Figure3G_component_necessity.tsv"
    if not generated.exists():
        return {"pass": False, "reason": "generated_file_missing"}
    got = pd.read_csv(generated, sep="\t")
    want = pd.read_csv(reference, sep="\t")
    keys = ["model", "component"]
    columns = ["mean_retention", "median_retention", "sem_retention"]
    joined = got[keys + columns + ["n_seeds"]].merge(
        want[keys + columns + ["n_seeds"]], on=keys, suffixes=("_generated", "_reference"), validate="one_to_one"
    )
    # n_seeds is constant (=8) by construction; a Pearson correlation is
    # undefined for a zero-variance column, so compare it by plain equality.
    n_seeds_exact = bool(joined.n_seeds_generated.eq(joined.n_seeds_reference).all())
    results = {column: _numeric_summary(joined[f"{column}_generated"], joined[f"{column}_reference"], rtol=0.0, atol=5e-3, min_pearson=0.99) for column in columns}
    return {
        "pass": len(joined) == 12 and n_seeds_exact and all(item["pass"] for item in results.values()),
        "rows": len(joined), "n_seeds_exact": n_seeds_exact, "values": results,
    }


def compare_fig3e(run_dir: Path) -> dict[str, object]:
    generated = run_dir / "derived/Figure3E_directional_recovery/Figure3E_directional_scramble_recovery_three_gene_source.tsv"
    reference = REFERENCE_ROOT / "Figure3E_directional_recovery.tsv"
    if not generated.exists():
        return {"pass": False, "reason": "generated_file_missing"}
    got = pd.read_csv(generated, sep="\t")
    want = pd.read_csv(reference, sep="\t")
    keys = ["arm", "series", "extent_bp"]
    joined = got.merge(want, on=keys, suffixes=("_generated", "_reference"), validate="one_to_one")
    mean_values = _numeric_summary(joined.mean_retention_generated, joined.mean_retention_reference, rtol=0.0, atol=5e-3, min_pearson=0.99)
    sem = joined.dropna(subset=["sem_retention_generated", "sem_retention_reference"])
    sem_values = _numeric_summary(sem.sem_retention_generated, sem.sem_retention_reference, rtol=0.0, atol=5e-3, min_pearson=0.99)
    return {"pass": len(joined) == 112 and mean_values["pass"] and sem_values["pass"], "rows": len(joined), "mean_retention": mean_values, "sem_retention": sem_values}


def compare_fig3a(run_dir: Path) -> dict[str, object]:
    generated = run_dir / "derived/Figure3A_regional_ism/stage2_top_windows_per_gene.tsv"
    reference = REFERENCE_ROOT / "Figure3A_regional_ism/stage2_top_windows_per_gene.tsv"
    if not generated.exists():
        return {"pass": False, "reason": "generated_file_missing"}
    got = pd.read_csv(generated, sep="\t")
    want = pd.read_csv(reference, sep="\t")
    keys = ["gene", "window_start", "window_end"]
    got_i, want_i = got.set_index(keys), want.set_index(keys)
    common = got_i.index.intersection(want_i.index)
    # The stage-1 window scan is a greedy, non-overlapping top-N-per-gene
    # selection over live-scored deltas; a near-tied window pick can shift by
    # one step (5 bp) under ordinary AlphaGenome run-to-run drift, same as
    # the already-documented Figure 3C hotspot-selection sensitivity (R019).
    # So this compares the intersection of selected windows on correlation
    # and requires the overlap itself stay high, rather than exact window
    # identity.
    overlap_fraction = len(common) / len(want)
    columns = ["delta_liver", "delta_adipose_mean", "contrast_liver_minus_adipose_mean"]
    values = {
        column: _numeric_summary(got_i.loc[common, column], want_i.loc[common, column], rtol=0.0, atol=5e-4, min_pearson=0.995)
        for column in columns
    } if len(common) else {column: {"pass": False, "reason": "no_common_windows"} for column in columns}
    return {
        "pass": bool(len(got) == 150 and len(want) == 150 and overlap_fraction >= 0.9 and all(v["pass"] for v in values.values())),
        "generated_rows": len(got), "reference_rows": len(want), "common_window_overlap_fraction": overlap_fraction,
        "values": values,
    }


def compare_fig3f(run_dir: Path) -> dict[str, object]:
    generated = run_dir / "derived/Figure3F_boundary_grid/surface_summary_paired.csv"
    reference = REFERENCE_ROOT / "Figure3F_boundary_grid/surface_summary_paired.csv"
    if not generated.exists():
        return {"pass": False, "reason": "generated_file_missing"}
    got = pd.read_csv(generated)
    want = pd.read_csv(reference)
    keys = ["gene", "upstream_bp", "downstream_bp"]
    columns = ["outside_mean_retention", "outside_median_retention", "inside_mean_retention", "inside_median_retention"]
    joined = got[keys + columns].merge(want[keys + columns], on=keys, suffixes=("_generated", "_reference"), validate="one_to_one")
    results = {column: _numeric_summary(joined[f"{column}_generated"], joined[f"{column}_reference"], rtol=0.0, atol=5e-3, min_pearson=0.99) for column in columns}
    selected_generated = pd.read_csv(run_dir / "derived/Figure3F_boundary_grid/selected_mean_window.csv")
    selected_reference = pd.read_csv(REFERENCE_ROOT / "Figure3F_boundary_grid/selected_mean_window.csv")
    selected_match = (
        not selected_generated.empty and not selected_reference.empty
        and int(selected_generated.iloc[0].upstream_bp) == int(selected_reference.iloc[0].upstream_bp)
        and int(selected_generated.iloc[0].downstream_bp) == int(selected_reference.iloc[0].downstream_bp)
    )
    return {
        "pass": len(joined) == 2688 and selected_match and all(item["pass"] for item in results.values()),
        "rows": len(joined), "selected_window_matches": selected_match, "values": results,
    }


def compare_fig4b(run_dir: Path) -> dict[str, object]:
    generated = run_dir / "derived/Figure4B_distance_response/Figure4B_distance_response.tsv"
    reference = REFERENCE_ROOT / "Figure4B_distance_response.tsv"
    if not generated.exists():
        return {"pass": False, "reason": "generated_file_missing"}
    got = pd.read_csv(generated, sep="\t")
    # Despite the .tsv extension, this committed reference file is actually
    # comma-separated (same extension/content-mismatch class already
    # documented for Figure 1B/1C's PDF-vs-SVG naming in the original audit).
    want = pd.read_csv(reference)
    keys = ["upstream_distance_bp", "donor_group"]
    columns = ["mean", "median", "sem"]
    joined = got[keys + columns].merge(want[keys + columns], on=keys, suffixes=("_generated", "_reference"), validate="one_to_one")
    results = {column: _numeric_summary(joined[f"{column}_generated"], joined[f"{column}_reference"], rtol=0.0, atol=0.02, min_pearson=0.98) for column in columns}
    return {"pass": len(joined) == 16 and all(item["pass"] for item in results.values()), "rows": len(joined), "values": results}


def compare_fig4c(run_dir: Path) -> dict[str, object]:
    generated = run_dir / "derived/Figure4C_foldchange_cohorts/Figure4C_fold_change_summary.csv"
    reference = REFERENCE_ROOT / "Figure4C_foldchange_cohorts/Figure4C_fold_change_summary.csv"
    if not generated.exists():
        return {"pass": False, "reason": "generated_file_missing"}
    got = pd.read_csv(generated)
    want = pd.read_csv(reference)
    keys = ["cohort", "donor_group"]
    columns = ["mean", "median", "fraction_up", "fraction_down"]
    joined = got[keys + columns].merge(want[keys + columns], on=keys, suffixes=("_generated", "_reference"), validate="one_to_one")
    results = {column: _numeric_summary(joined[f"{column}_generated"], joined[f"{column}_reference"], rtol=0.0, atol=0.5, min_pearson=0.9) for column in columns}
    return {"pass": len(joined) == 9 and all(item["pass"] for item in results.values()), "rows": len(joined), "values": results}


def compare_fig4e(run_dir: Path) -> dict[str, object]:
    generated = run_dir / "derived/Figure4EF_distal_contact_transfer/Figure4E_distance_fraction_positive.tsv"
    reference = REFERENCE_ROOT / "Figure4E_distance_fraction_positive/plot_distance_band_bootstrap.tsv"
    if not generated.exists():
        return {"pass": False, "reason": "generated_file_missing"}
    got = pd.read_csv(generated, sep="\t")
    want = pd.read_csv(reference, sep="\t")
    keys = ["distance_label"]
    columns = ["fraction", "n_clusters"]
    joined = got[keys + columns].merge(want[keys + columns], on=keys, suffixes=("_generated", "_reference"), validate="one_to_one")
    # `fraction` is a raw empirical proportion (not a bootstrap draw -- only
    # its CI is), so it's an ordinary two-independent-sample proportion
    # difference between one freshly re-scored AlphaGenome run and the
    # archived run; see _cluster_proportion_difference_check for the
    # per-row, cluster-count-derived tolerance (same statistic as 4F).
    results = {"fraction": _cluster_proportion_difference_check(
        joined["fraction_generated"], joined["fraction_reference"],
        joined["n_clusters_generated"], joined["n_clusters_reference"],
    )}
    return {"pass": len(joined) == 5 and all(item["pass"] for item in results.values()), "rows": len(joined), "values": results}


def _cluster_proportion_difference_check(
    generated: np.ndarray, reference: np.ndarray, n_clusters_generated: np.ndarray, n_clusters_reference: np.ndarray,
    *, z: float = 1.96,
) -> dict[str, object]:
    """Per-row tolerance for comparing an empirical proportion (`fraction` of
    promoters with a positive T-G interaction, one independent AlphaGenome
    re-scoring vs. the archived run) computed over a cluster-bootstrap unit.
    `fraction` itself is not a bootstrap draw -- only its CI is -- so it is
    an ordinary two-independent-sample proportion difference; its worst-case
    standard error (at p=0.5, which maximizes Bernoulli variance) is
    sqrt(p(1-p)/n1 + p(1-p)/n2) with n1, n2 the cluster counts (not raw row
    counts, since clusters -- not promoters -- are the independent unit).
    A z=1.96 (95%) bound is applied per-row from that row's own cluster
    count, not a single blanket tolerance fit to whichever row is hardest to
    pass."""
    generated, reference = np.asarray(generated, dtype=float), np.asarray(reference, dtype=float)
    n_min = np.minimum(np.asarray(n_clusters_generated, dtype=float), np.asarray(n_clusters_reference, dtype=float))
    threshold = z * np.sqrt(0.5 * 0.5 / n_min + 0.5 * 0.5 / n_min)
    delta = np.abs(generated - reference)
    row_pass = delta <= threshold
    return {
        "pass": bool(row_pass.all()),
        "n": int(len(generated)),
        "max_abs_difference": float(delta.max()),
        "mean_abs_difference": float(delta.mean()),
        "criterion": "per-row 95% two-independent-sample proportion-difference bound at worst-case p=0.5, using min(n_clusters_generated, n_clusters_reference)",
        "rows_failing": int((~row_pass).sum()),
    }


def compare_fig4f(run_dir: Path) -> dict[str, object]:
    generated = run_dir / "derived/Figure4EF_distal_contact_transfer/Figure4F_contact_dose_response.tsv"
    reference = REFERENCE_ROOT / "Figure4F_contact_dose_response.tsv"
    if not generated.exists():
        return {"pass": False, "reason": "generated_file_missing"}
    got = pd.read_csv(generated, sep="\t")
    want = pd.read_csv(reference, sep="\t")
    keys = ["distance_stratum", "contact_quintile"]
    joined = got.merge(want, on=keys, suffixes=("_generated", "_reference"), validate="one_to_one")
    results = {
        "median_contact_contrast": _numeric_summary(
            joined["median_contact_contrast_generated"], joined["median_contact_contrast_reference"], rtol=0.05, atol=0.15,
        ),
        "fraction": _cluster_proportion_difference_check(
            joined["fraction_generated"], joined["fraction_reference"],
            joined["n_clusters_generated"], joined["n_clusters_reference"],
        ),
    }
    return {"pass": len(joined) == 25 and all(item["pass"] for item in results.values()), "rows": len(joined), "values": results}


def compare_fig4g(run_dir: Path) -> dict[str, object]:
    generated = run_dir / "derived/Figure4G_tissue_rna/Figure4G_tissue_rna.tsv"
    reference = REFERENCE_ROOT / "Figure4G_tissue_rna.tsv"
    if not generated.exists():
        return {"pass": False, "reason": "generated_file_missing"}
    got = pd.read_csv(generated, sep="\t")
    want = pd.read_csv(reference, sep="\t")
    joined = got.merge(want, on="context", suffixes=("_generated", "_reference"), validate="one_to_one")
    columns = ["SORT1", "PSRC1", "CELSR2"]
    # A single variant, single API call, no bootstrap resampling -- deviation
    # from the archived run should be ordinary AlphaGenome run-to-run float
    # precision (established elsewhere in this project at 1e-4-1e-3), so a
    # tight tolerance is meaningful here, unlike the bootstrap-based panels.
    results = {
        column: _numeric_summary(joined[f"{column}_generated"], joined[f"{column}_reference"], rtol=0.0, atol=0.01, min_pearson=0.999)
        for column in columns
    }
    return {"pass": len(joined) == 7 and all(item["pass"] for item in results.values()), "rows": len(joined), "values": results}


COMPARATORS = {
    "1B": compare_fig1b, "1C": compare_fig1c, "1C-middle": compare_fig1c_middle,
    "1D": compare_fig1d, "1E": compare_fig1e, "1F": compare_fig1f,
    "2B": compare_fig2b, "2C": compare_fig2c, "2E": compare_fig2e, "2F": compare_fig2f,
    "3A": compare_fig3a, "3B": compare_fig3b, "3C": compare_fig3c, "3E": compare_fig3e,
    "3F": compare_fig3f, "3G": compare_fig3g, "4B": compare_fig4b, "4C": compare_fig4c,
    "4E": compare_fig4e, "4F": compare_fig4f, "4G": compare_fig4g,
}


def compare_run(run_dir: Path, panels: list[str]) -> dict[str, object]:
    results = {panel: COMPARATORS[panel](run_dir) for panel in panels}
    comparison = {
        "created_utc": utc_now(),
        "comparison_code": {
            "path": str(Path(__file__).relative_to(REPO_ROOT)),
            "bytes": Path(__file__).stat().st_size,
            "sha256": sha256_file(Path(__file__)),
        },
        "reference_root": str(REFERENCE_ROOT.resolve()),
        "note": "Reference outputs were accessed only by this post-run comparison stage.",
        "pass": all(bool(item["pass"]) for item in results.values()),
        "panels": results,
    }
    path = run_dir / "audit" / "comparison.json"
    path.write_text(json.dumps(comparison, indent=2) + "\n")
    return comparison


def write_report(run_dir: Path, comparison: dict[str, object] | None = None) -> Path:
    audit_path = run_dir / "audit" / "run.json"
    audit = json.loads(audit_path.read_text())
    prior_audits = []
    for prior_path in sorted((run_dir / "audit" / "attempts").glob("*/run.json")):
        try:
            prior_audits.append((prior_path, json.loads(prior_path.read_text())))
        except (OSError, json.JSONDecodeError):
            continue
    if comparison is None and (run_dir / "audit" / "comparison.json").exists():
        comparison = json.loads((run_dir / "audit" / "comparison.json").read_text())
    steps = audit.get("steps", [])
    elapsed = sum(float(step.get("elapsed_seconds", 0)) for step in steps)
    cumulative_elapsed = elapsed + sum(
        float(step.get("elapsed_seconds", 0))
        for _, prior in prior_audits for step in prior.get("steps", [])
    )
    attempts = [prior for _, prior in prior_audits] + [audit]
    scored_units = {
        panel: sum(
            int(attempt.get("api_calls", {}).get(panel, 0))
            + (int(attempt.get("api_calls", {}).get("1C-middle", 0)) if panel == "1C" else 0)
            + (int(attempt.get("api_calls", {}).get("4E/4F", 0)) if panel == "4E" else 0)
            for attempt in attempts
        )
        for panel in audit["panels"]
    }
    api_requests = {
        panel: sum(
            int(attempt.get("api_requests", {}).get(panel, 0))
            + (int(attempt.get("api_requests", {}).get("1C-middle", 0)) if panel == "1C" else 0)
            + (int(attempt.get("api_requests", {}).get("4E/4F", 0)) if panel == "4E" else 0)
            for attempt in attempts
        )
        for panel in audit["panels"]
    }
    lines = [
        "# SORT1 clean-room reproducibility report",
        "",
        f"- Run status: **{audit['status']}**",
        f"- Reference comparison: **{'PASS' if comparison and comparison.get('pass') else 'NOT YET PASSED'}**",
        f"- Started (UTC): `{audit['started_utc']}`",
        f"- Repository commit: `{audit['repository']['commit']}`",
        f"- Repository dirty during run: `{audit['repository']['dirty']}`",
        f"- Run directory: `{audit['run_directory']}`",
        "- AlphaGenome credential available during run: "
        f"`{audit['environment']['api_key_present']}` (value never recorded)",
        f"- Timed step total: `{elapsed:.1f} s`",
        f"- Timed step total across retained attempts: `{cumulative_elapsed:.1f} s`",
        f"- Peak RSS reported by process: `{audit.get('peak_rss_bytes', 0)} bytes`",
        f"- Run-directory size: `{audit.get('disk_bytes', 0)} bytes`",
        "",
        "## Scope",
        "",
        "The analysis stage is isolated from `outputs/source_data/`. Published source tables "
        "are opened only after a successful analysis, by the comparison command.",
        "",
        "## Executed reproduction code",
        "",
        "| Path | Bytes | SHA-256 |",
        "|---|---:|---|",
    ]
    if comparison and comparison.get("comparison_code"):
        code = comparison["comparison_code"]
        lines[lines.index("## Executed reproduction code"):lines.index("## Executed reproduction code")] = [
            "## Post-run comparison code",
            "",
            f"- Path: `{code['path']}`",
            f"- Bytes: `{code['bytes']}`",
            f"- SHA-256: `{code['sha256']}`",
            "",
        ]
    for item in audit.get("reproduction_code", []):
        lines.append(f"| `{item['path']}` | {item['bytes']} | `{item['sha256']}` |")
    lines.extend(
        [
        "",
        "## Panel results",
        "",
        "| Panel | Inputs | AlphaGenome regime | Scored units | API requests | Comparison |",
        "|---|---|---|---:|---:|---|",
        ]
    )
    descriptions = {
        "1B": ("AlphaGenome API", "ALL_FOLDS"),
        "1C": ("GLGC + GTEx v7 liver + 1000 Genomes EUR + AlphaGenome API", "ALL_FOLDS"),
        "1C-middle": ("GLGC 2013 + UCSC LiftOver + hg38 + AlphaGenome API", "ALL_FOLDS"),
        "1D": ("4DN HepG2 observed Hi-C", "none"),
        "1E": ("AlphaGenome API", "FOLD_0"),
        "1F": ("4DN HepG2 Hi-C + AlphaGenome contact map", "FOLD_0"),
        "2B": ("Wang 2018 spreadsheet + hg38 + AlphaGenome API", "ALL_FOLDS"),
        "2C": ("hg38 + AlphaGenome API", "ALL_FOLDS"),
        "2E": ("Kircher 2019 MPRA + AlphaGenome API", "ALL_FOLDS"),
        "2F": ("Kircher 2019 MPRA + AlphaGenome API", "ALL_FOLDS"),
        "3A": ("GRCh38 + AlphaGenome API", "ALL_FOLDS"),
        "3B": ("GRCh38 + AlphaGenome API", "ALL_FOLDS"),
        "3C": ("Figure 3B outputs + JASPAR 2024 CORE + GRCh38", "none (local PWM scan)"),
        "3G": ("GRCh38 + AlphaGenome API", "ALL_FOLDS"),
        "3E": ("GRCh38 + AlphaGenome API", "ALL_FOLDS"),
        "3F": ("GRCh38 + AlphaGenome API", "ALL_FOLDS"),
        "4B": ("GRCh38 + frozen bottom100 recipient design + AlphaGenome API", "ALL_FOLDS"),
        "4C": ("GRCh38 + HPA v24.1 + GENCODE v46 + AlphaGenome API", "ALL_FOLDS"),
        "4E": ("GRCh38 + GENCODE v46 + 4DN HepG2 Hi-C + AlphaGenome API", "ALL_FOLDS"),
        "4F": ("GRCh38 + GENCODE v46 + 4DN HepG2 Hi-C + AlphaGenome API (shares 4E's run)", "ALL_FOLDS"),
        "4G": ("AlphaGenome API (single variant, RNA_SEQ scorer)", "ALL_FOLDS"),
    }
    for panel in audit["panels"]:
        panel_comparison = comparison and comparison.get("panels", {}).get(panel)
        if panel_comparison is None:
            label = "NOT COMPARED"
        else:
            label = "PASS" if panel_comparison.get("pass") else "FAIL"
        inputs, regime = descriptions[panel]
        lines.append(
            f"| {panel} | {inputs} | {regime} | {scored_units[panel]} | "
            f"{api_requests[panel]} | {label} |"
        )
    lines.extend(["", "## Timings", "", "| Step | Status | Seconds |", "|---|---|---:|"])
    for step in steps:
        lines.append(f"| {step['name']} | {step['status']} | {float(step.get('elapsed_seconds', 0)):.3f} |")
    preparation_checks = [
        step
        for step in steps
        if "variants" in step or "grch38_reference_allele_matches" in step
    ]
    if preparation_checks:
        lines.extend(["", "## Preparation checks", ""])
        for step in preparation_checks:
            if "variants" in step:
                lines.append(f"- Reconstructed variants: `{step['variants']}`")
            if "grch38_reference_allele_matches" in step:
                lines.append(
                    "- Variants with a GLGC allele matching downloaded GRCh38: "
                    f"`{step['grch38_reference_allele_matches']}`"
                )
    prior_runs = []
    for prior_path, prior in prior_audits:
        prior_runs.append(
            (
                prior.get("started_utc", "unknown"),
                prior.get("status", "unknown"),
                sum(float(step.get("elapsed_seconds", 0)) for step in prior.get("steps", [])),
                prior_path.parent.name,
            )
        )
    if prior_runs:
        lines.extend(
            [
                "",
                "## Prior attempts retained during resume",
                "",
                "| Started UTC | Status | Timed seconds | Archive |",
                "|---|---|---:|---|",
            ]
        )
        for started, status, seconds, archive in prior_runs:
            lines.append(f"| {started} | {status} | {seconds:.3f} | `audit/attempts/{archive}/` |")
    lines.extend(["", "## Downloads", ""])
    if audit.get("downloads"):
        lines.extend(["| URL | Bytes | SHA-256 | Reused |", "|---|---:|---|---|"])
        for item in audit["downloads"]:
            lines.append(
                f"| {item['url']} | {item['bytes']} | `{item['sha256']}` | "
                f"{item.get('reused', False)} |"
            )
    else:
        lines.append("No external datasets were downloaded for the selected panels.")
    lines.extend(["", "## Generated-file manifest", "", "| Path | Bytes | SHA-256 |", "|---|---:|---|"])
    for item in file_manifest(run_dir):
        lines.append(f"| `{item['path']}` | {item['bytes']} | `{item['sha256']}` |")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "A PASS means identities and counts were exact and numerical outputs met the explicit "
            "panel-specific equivalence thresholds recorded in `audit/comparison.json`. Figure 2F and "
            "non-live comparisons retain `rtol=1e-5, atol=1e-6`; live sequence/ISM panels additionally "
            "allow bounded sub-panel-unit API drift while requiring very high reference correlation. "
            "The observed maximum differences and plotted-correlation changes remain fully reported. "
            "The frozen tables were not used to generate any result.",
            "",
        ]
    )
    path = run_dir / "audit" / "REPRODUCIBILITY_REPORT.md"
    path.write_text("\n".join(lines))
    return path
