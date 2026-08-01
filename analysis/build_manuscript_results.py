#!/usr/bin/env python3
"""Rebuild the headline manuscript-result ledger from release source tables."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "outputs" / "source_data"
OUT = ROOT / "outputs" / "manuscript_results.tsv"


def main() -> None:
    rows: list[dict[str, object]] = []

    def add(result_id: str, analysis_id: str, metric: str, value: object, unit: str,
            model: str, population_or_n: str, source_table: str, notes: str = "") -> None:
        rows.append({
            "result_id": result_id,
            "analysis_id": analysis_id,
            "metric": metric,
            "value": value,
            "unit": unit,
            "model": model,
            "population_or_n": population_or_n,
            "source_table": source_table,
            "status": "verified_from_release_table",
            "notes": notes,
        })

    # Figure 1: direct and LD-tagging effects.
    fig1_path = "outputs/source_data/Figure1C_eqtl_direct_tagging.tsv"
    fig1 = pd.read_csv(ROOT / fig1_path, sep="\t")
    rs = fig1.loc[fig1["rsid"].eq("rs12740374")].iloc[0]
    for gene in ("SORT1", "CELSR2", "PSRC1"):
        direct_value = rs[f"ag_single_snp_term_filled_{gene}"]
        add(f"fig1c_rs127_direct_{gene.lower()}", "fig1_eqtl_prioritization",
            f"rs12740374 direct RNA effect, {gene}", direct_value,
            "GeneMaskLFCScorer value", "ALL_FOLDS", "rs12740374", fig1_path)
        add(f"fig1c_rs127_percent_{gene.lower()}", "fig1_eqtl_prioritization",
            f"rs12740374 implied RNA percentage change, {gene}", 100 * np.expm1(direct_value),
            "percent", "ALL_FOLDS", "rs12740374", fig1_path,
            "Calculated as 100 * expm1(GeneMaskLFCScorer value).")
        observed = fig1[f"eqtl_liver_{gene}"]
        direct = fig1[f"ag_single_snp_term_filled_{gene}"]
        tagging = fig1[f"ag_model_snp_plus_covvar_rs127_for_plot_{gene}"]
        valid = observed.notna() & direct.notna()
        add(f"fig1c_direct_pearson_{gene.lower()}", "fig1_eqtl_prioritization",
            f"Observed eQTL versus direct score Pearson r, {gene}", pearsonr(observed[valid], direct[valid]).statistic,
            "Pearson r", "ALL_FOLDS", f"n={int(valid.sum())} variants", fig1_path)
        valid = observed.notna() & tagging.notna()
        add(f"fig1c_tagging_pearson_{gene.lower()}", "fig1_eqtl_prioritization",
            f"Observed eQTL versus rs12740374-tagging model Pearson r, {gene}",
            pearsonr(observed[valid], tagging[valid]).statistic, "Pearson r", "ALL_FOLDS",
            f"n={int(valid.sum())} variants", fig1_path,
            "Tagging model explicitly assumes rs12740374 is the driver.")
        add(f"fig1c_rs127_absolute_rank_{gene.lower()}", "fig1_eqtl_prioritization",
            f"rs12740374 absolute direct-effect rank, {gene}",
            int(direct.abs().rank(method="min", ascending=False).loc[fig1["rsid"].eq("rs12740374")].iloc[0]),
            "rank", "ALL_FOLDS", f"n={int(direct.notna().sum())} variants", fig1_path)

    contact_path = "outputs/source_data/Figure1F_promoter_contact_percentiles.tsv"
    contact = pd.read_csv(ROOT / contact_path, sep="\t")
    for _, row in contact.iterrows():
        prefix = "observed" if row["scale"].startswith("Observed") else "fold0"
        gene = row["gene"]
        add(f"fig1f_{prefix}_contact_{gene.lower()}", "fig1_contact_percentiles",
            f"{row['scale']} contact value, {gene}", row["contact_value"],
            "observed/expected or model contact value", "EXPERIMENTAL" if prefix == "observed" else "FOLD_0",
            gene, contact_path)
        add(f"fig1f_{prefix}_percentile_{gene.lower()}", "fig1_contact_percentiles",
            f"{row['scale']} same-distance percentile, {gene}", row["same_distance_percentile"],
            "percentile", "EXPERIMENTAL" if prefix == "observed" else "FOLD_0", gene, contact_path)

    # Figure 2: repair outcomes and Kircher benchmarks.
    repair_path = "outputs/source_data/Figure2B_top50_repair_outcomes/matrix.csv"
    repair = pd.read_csv(ROOT / repair_path).set_index("gene_symbol")
    for gene, values in repair.iterrows():
        add(f"fig2b_negative_count_{gene.lower()}", "fig2_repair_top50",
            f"Repair products with negative predicted RNA effect, {gene}", int((values < 0).sum()),
            "of 50 products", "ALL_FOLDS", "50 ranked human/mouse products", repair_path)
        add(f"fig2b_median_{gene.lower()}", "fig2_repair_top50",
            f"Median predicted RNA change, {gene}", values.median(), "percent", "ALL_FOLDS",
            "50 ranked human/mouse products", repair_path)

    kircher_path = "outputs/source_data/Figure2E_kircher_correlations.tsv"
    kircher = pd.read_csv(ROOT / kircher_path, sep="\t")
    for label, column in {
        "three_gene_mean_RNA": "ag_rna_3gene_mean_percent_change",
        "SORT1_RNA": "ag_rna_percent_change_SORT1",
        "PSRC1_RNA": "ag_rna_percent_change_PSRC1",
        "CELSR2_RNA": "ag_rna_percent_change_CELSR2",
        "ATAC": "ag_atac_mean_score",
        "H3K27ac": "ag_h3k27ac_mean_score",
    }.items():
        valid = kircher["kircher_primary_log2_effect"].notna() & kircher[column].notna()
        add(f"fig2e_pearson_{label.lower()}", "fig2_kircher_sort1",
            f"Kircher MPRA versus AlphaGenome Pearson r, {label}",
            pearsonr(kircher.loc[valid, "kircher_primary_log2_effect"], kircher.loc[valid, column]).statistic,
            "Pearson r", "ALL_FOLDS", f"n={int(valid.sum())} substitutions", kircher_path)
        add(f"fig2e_spearman_{label.lower()}", "fig2_kircher_sort1",
            f"Kircher MPRA versus AlphaGenome Spearman rho, {label}",
            spearmanr(kircher.loc[valid, "kircher_primary_log2_effect"], kircher.loc[valid, column]).statistic,
            "Spearman rho", "ALL_FOLDS", f"n={int(valid.sum())} substitutions", kircher_path)

    multi_path = "outputs/source_data/Figure2F_kircher_multielement/RESULTS.json"
    multi = json.loads((ROOT / multi_path).read_text())
    for modality in ("ATAC", "DNASE"):
        result = multi["random_effects_by_modality"][modality]
        add(f"fig2f_pooled_{modality.lower()}", "fig2_kircher_multielement",
            f"Random-effects pooled Spearman rho, {modality}", result["pooled_spearman_rho"],
            "Spearman rho", "ALL_FOLDS", f"k={result['k_elements']} elements", multi_path,
            f"95% CI {result['pooled_95ci_low']:.6g} to {result['pooled_95ci_high']:.6g}; I2={result['I2_percent']:.3f}%")
        for metric, unit in (("Q", "Cochran Q"), ("Q_p", "P value"), ("I2_percent", "percent")):
            add(f"fig2f_{modality.lower()}_{metric.lower()}", "fig2_kircher_multielement",
                f"Random-effects heterogeneity {metric}, {modality}", result[metric], unit,
                "ALL_FOLDS", f"k={result['k_elements']} elements", multi_path)
    for modality, result in multi["matched_four_random_effects_by_modality"].items():
        add(f"fig2f_matched4_{modality.lower()}", "fig2_kircher_multielement",
            f"Matched-four random-effects Spearman rho, {modality}", result["pooled_spearman_rho"],
            "Spearman rho", "ALL_FOLDS", "k=4 matched elements", multi_path,
            f"95% CI {result['pooled_95ci_low']:.6g} to {result['pooled_95ci_high']:.6g}")
    element_path = "outputs/source_data/Figure2F_kircher_multielement/kircher_multielement_element_statistics.tsv"
    elements = pd.read_csv(ROOT / element_path, sep="\t")
    for _, row in elements.iterrows():
        add(f"fig2f_{row['element'].lower()}_{row['modality'].lower()}_rho", "fig2_kircher_multielement",
            f"Element-level Spearman rho, {row['element']}, {row['modality']}", row["spearman_rho"],
            "Spearman rho", "ALL_FOLDS", f"n={int(row['n'])} substitutions", element_path,
            f"Cell context: {row['cell_line']}; 95% CI {row['spearman_ci_low']:.6g} to {row['spearman_ci_high']:.6g}")

    # Figure 3: frozen boundary and component results.
    for criterion, filename in (("mean", "selected_mean_window.csv"), ("median", "selected_median_window.csv")):
        path = f"outputs/source_data/Figure3F_boundary_grid/{filename}"
        row = pd.read_csv(ROOT / path).iloc[0]
        for metric, unit in (("window_length_bp", "bp"), ("outside_mean_retention", "fraction"),
                             ("outside_median_retention", "fraction"), ("inside_mean_retention", "fraction")):
            add(f"fig3f_{criterion}_{metric}", "fig3_boundary_grid",
                f"Shortest {criterion}-passing window {metric}", row[metric], unit, "ALL_FOLDS",
                "8 frozen scramble seeds", path)

    component_path = "outputs/source_data/Figure3G_component_necessity.tsv"
    component = pd.read_csv(ROOT / component_path, sep="\t")
    for _, row in component.iterrows():
        add(f"fig3g_retention_{row['component'].lower().replace('/', '_').replace(' ', '_')}",
            "fig3_component_necessity", f"Mean three-gene retention after scrambling {row['component']}",
            row["mean_retention"], "fraction", row["model"], f"n={int(row['n_seeds'])} seeds", component_path)

    recovery_path = "outputs/source_data/Figure3E_directional_recovery.tsv"
    recovery = pd.read_csv(ROOT / recovery_path, sep="\t")
    for arm in ("Upstream", "Downstream"):
        part = recovery.loc[(recovery["series"] == "3-gene mean") & (recovery["arm"] == arm)]
        for label, row in (("protected_only", part.loc[part["extent_bp"].eq(0)].iloc[0]),
                           ("maximum", part.loc[part["mean_retention"].idxmax()])):
            add(f"fig3e_{arm.lower()}_{label}", "fig3_directional_recovery",
                f"Three-gene mean retention, {arm.lower()} restoration, {label}", row["mean_retention"],
                "fraction", "ALL_FOLDS", "8 frozen scramble seeds", recovery_path,
                f"Native sequence restored: {int(row['extent_bp'])} bp")

    audit_path = "outputs/source_data/FigureS7_native_sequence_audit.tsv"
    audit = pd.read_csv(ROOT / audit_path, sep="\t")
    audit = audit.loc[audit["construct"].ne("intact_T")]
    paired = audit.pivot_table(index=["construct", "gene"], columns="model", values="loss_vs_intact_T").dropna()
    add("figs7_direction_agreement_all_genes", "figs7_hotspot_heldout_audit",
        "ALL_FOLDS/FOLD_0 direction agreement across gene-window comparisons",
        int((paired["ALL_FOLDS"] * paired["FOLD_0"] > 0).sum()), "of 21 comparisons",
        "ALL_FOLDS;FOLD_0", "7 windows x 3 genes", audit_path)
    sort1_audit = paired.reset_index().loc[lambda x: x["gene"].eq("SORT1")]
    add("figs7_sort1_direction_agreement", "figs7_hotspot_heldout_audit",
        "ALL_FOLDS/FOLD_0 direction agreement across SORT1 windows",
        int((sort1_audit["ALL_FOLDS"] * sort1_audit["FOLD_0"] > 0).sum()), "of 7 windows",
        "ALL_FOLDS;FOLD_0", "7 frozen windows", audit_path)
    add("figs7_sort1_rank_concordance", "figs7_hotspot_heldout_audit",
        "ALL_FOLDS/FOLD_0 SORT1 window-rank concordance",
        spearmanr(sort1_audit["ALL_FOLDS"], sort1_audit["FOLD_0"]).statistic,
        "Spearman rho", "ALL_FOLDS;FOLD_0", "7 frozen windows", audit_path)

    # Figure 4: promoter transfer and distal contact analysis.
    fig4b_path = "outputs/source_data/Figure4B_distance_response.tsv"
    fig4b = pd.read_csv(ROOT / fig4b_path, sep=None, engine="python")
    for _, row in fig4b.loc[fig4b["upstream_distance_bp"].eq(30)].iterrows():
        add(f"fig4b_30bp_{row['donor_group']}", "fig4_transfer_distance",
            f"Mean RNA change at 30 bp, {row['donor_group']}", row["mean"], "AlphaGenome RNA delta",
            "ALL_FOLDS", f"n={int(row['n'])} recipients", fig4b_path)

    fig4c_path = "outputs/source_data/Figure4C_foldchange_cohorts/Figure4C_fold_change_summary.csv"
    fig4c = pd.read_csv(ROOT / fig4c_path)
    for _, row in fig4c.iterrows():
        key = f"{row['cohort']}_{row['donor_group']}"
        add(f"fig4c_median_{key}", "fig4_foldchange_cohorts",
            f"Median log2 construct/native effect, {row['cohort']}, {row['donor_group']}", row["median"],
            "log2 fold change", "ALL_FOLDS", f"n={int(row['n'])} recipients", fig4c_path)
        add(f"fig4c_fraction_up_{key}", "fig4_foldchange_cohorts",
            f"Fraction above native, {row['cohort']}, {row['donor_group']}", row["fraction_up"],
            "fraction", "ALL_FOLDS", f"n={int(row['n'])} recipients", fig4c_path)

    fig4e_path = "outputs/source_data/Figure4E_distance_fraction_positive/plot_distance_band_bootstrap.tsv"
    fig4e = pd.read_csv(ROOT / fig4e_path, sep="\t")
    for _, row in fig4e.iterrows():
        add(f"fig4e_positive_{row['distance_label'].replace('–', '_').replace(' ', '')}",
            "fig4_distance_fraction_positive",
            f"High-contact site stronger fraction, {row['distance_label']}", row["fraction"],
            "fraction", "ALL_FOLDS", f"n={int(row['n'])} promoters", fig4e_path,
            f"95% bootstrap CI {row['ci_low']:.6g} to {row['ci_high']:.6g}")

    fig4f_path = "outputs/source_data/Figure4F_contact_dose_response.tsv"
    fig4f = pd.read_csv(ROOT / fig4f_path, sep="\t")
    for _, row in fig4f.iterrows():
        key = f"{row['distance_stratum']}_q{int(row['contact_quintile'])}".replace("–", "_").replace(",", "")
        add(f"fig4f_fraction_{key}", "fig4_contact_dose_response",
            f"High-contact site stronger fraction, {row['distance_stratum']}, contact quintile {int(row['contact_quintile'])}",
            row["fraction"], "fraction", "ALL_FOLDS", f"n={int(row['n'])} promoters", fig4f_path,
            f"Median observed contact contrast {row['median_contact_contrast']:.6g}; 95% bootstrap CI {row['ci_low']:.6g} to {row['ci_high']:.6g}")

    tissue_path = "outputs/source_data/Figure4G_tissue_rna.tsv"
    tissue = pd.read_csv(ROOT / tissue_path, sep="\t").set_index("context")
    for gene in ("SORT1", "PSRC1", "CELSR2"):
        add(f"fig4g_liver_{gene.lower()}", "fig4_tissue_rna",
            f"rs12740374 liver RNA delta, {gene}", tissue.loc["liver", gene],
            "AlphaGenome RNA delta", "ALL_FOLDS", "liver ontology", tissue_path)

    # Figure S1: contact summarization, allele effects, and cross-context baseline architecture.
    s1b_path = "outputs/source_data/FigureS1B_tss_bin_contact.tsv"
    s1b = pd.read_csv(ROOT / s1b_path, sep="\t")
    for _, row in s1b.iterrows():
        regime = "EXPERIMENTAL" if str(row["source"]).startswith("Experimental") else "FOLD_0"
        add(f"figs1b_{regime.lower()}_{row['gene'].lower()}_percentile", "figs1_tss_bin_contact",
            f"TSS-bin same-distance contact percentile, {row['source']}, {row['gene']}",
            row["same_distance_percentile"], "percentile", regime,
            f"same-distance null n={int(row['same_distance_null_n'])}", s1b_path)
    s1c_path = "outputs/source_data/FigureS1C_contact_allele_delta/panelB_promoter_allele_delta.tsv"
    s1c = pd.read_csv(ROOT / s1c_path, sep="\t")
    for _, row in s1c.iterrows():
        add(f"figs1c_percent_change_{row['gene'].lower()}", "figs1_contact_allele_delta",
            f"ALT versus REF promoter-contact change, {row['gene']}", row["percent_change_implied_observed_expected"],
            "percent", "ALL_FOLDS", row["gene"], s1c_path)
    s1d_path = "outputs/source_data/FigureS1D_contact_contexts.tsv"
    s1d = pd.read_csv(ROOT / s1d_path, sep="\t")
    context_pivot = s1d.pivot(index="ontology_curie", columns="gene", values="same_distance_percentile")
    sort1_highest = (context_pivot["SORT1"] > context_pivot[["CELSR2", "PSRC1"]].max(axis=1)).sum()
    add("figs1d_sort1_highest_contexts", "figs1_contact_contexts",
        "Contexts in which SORT1 had the highest same-distance contact percentile", int(sort1_highest),
        f"of {len(context_pivot)} ontologies", "ALL_FOLDS", f"n={len(context_pivot)} ontologies", s1d_path)

    # Figure S2: direct-score stability when the locus is held out.
    s2_path = "outputs/source_data/FigureS2_eqtl_fold0/figureS1C_all_folds_vs_fold0_source.tsv"
    s2 = pd.read_csv(ROOT / s2_path, sep="\t")
    for gene, part in s2.groupby("gene"):
        add(f"figs2_rank_concordance_{gene.lower()}", "figs2_eqtl_fold0",
            f"ALL_FOLDS versus FOLD_0 direct-score rank concordance, {gene}",
            spearmanr(part["ALL_FOLDS"], part["FOLD_0"]).statistic, "Spearman rho",
            "ALL_FOLDS;FOLD_0", f"n={len(part)} variants", s2_path)
        rsrow = part.loc[part["rsid"].eq("rs12740374")].iloc[0]
        for regime in ("ALL_FOLDS", "FOLD_0"):
            add(f"figs2_rs127_{regime.lower()}_{gene.lower()}", "figs2_eqtl_fold0",
                f"rs12740374 direct score, {gene}, {regime}", rsrow[regime],
                "GeneMaskLFCScorer value", regime, "rs12740374", s2_path)

    # Figure S3: caQTL direct and tagging-model correlations.
    s3_path = "outputs/source_data/FigureS3_caqtl_locus/figureS2_correlations.tsv"
    s3 = pd.read_csv(ROOT / s3_path, sep="\t")
    for _, row in s3.iterrows():
        for metric in ("pearson_r", "spearman_rho"):
            add(f"figs3_{row['peak_group']}_{row['regime'].lower()}_{row['model']}_{metric}",
                "figs3_caqtl_locus", f"Observed caQTL versus {row['model']} model {metric}, {row['peak_group']}, {row['regime']}",
                row[metric], "Pearson r" if metric == "pearson_r" else "Spearman rho", row["regime"],
                f"n={int(row['n'])} variants", s3_path,
                "Tagging results assume rs12740374 is the driver and are explanatory, not an independent causal test." if row["model"] == "tagging" else "")

    # Figure S4: human repair depth and conditional frequency-weighted mixture.
    s4a_path = "outputs/source_data/FigureS4A_sequencing_depth.tsv"
    s4a = pd.read_csv(ROOT / s4a_path, sep="\t").iloc[0]
    for metric, unit in (("total_reads", "reads"), ("modeled_simple_indel_reads", "reads"),
                         ("editing_fraction_percent", "percent"), ("modeled_indel_coverage_percent", "percent")):
        add(f"figs4a_{metric}", "figs4_human_repair_depth", metric.replace("_", " "), s4a[metric], unit,
            "EXPERIMENTAL", "primary human hepatocytes", s4a_path)
    s4b_path = "outputs/source_data/FigureS4B_freqweighted_heatmap/FigureS4B_frequency_weighted_all_repair_outcomes_summary.csv"
    s4b = pd.read_csv(ROOT / s4b_path)
    for _, row in s4b.iterrows():
        add(f"figs4b_weighted_{row['gene'].lower()}", "figs4_human_repair_weighted",
            f"Conditional frequency-weighted RNA change, {row['gene']}", row["frequency_weighted_mixture_percent"],
            "percent", "ALL_FOLDS", f"n={int(row['n_human_products'])} human repair products", s4b_path,
            "Conditioned on the 1,396 modeled simple-indel reads; not multiplied by the bulk editing fraction.")

    # Figure S5: default-versus-held-out SORT1-locus MPRA benchmarks.
    s5_specs = {
        "RNA_substitutions": ("FigureS5A_RNA_substitutions.tsv", {
            "three_gene_mean": "ag_rna_3gene_mean_percent_change", "SORT1": "ag_rna_percent_change_SORT1",
            "PSRC1": "ag_rna_percent_change_PSRC1", "CELSR2": "ag_rna_percent_change_CELSR2"}),
        "ATAC_substitutions": ("FigureS5B_ATAC_substitutions.tsv", {"ATAC": "ag_atac_mean_score"}),
        "H3K27ac_substitutions": ("FigureS5C_H3K27ac_substitutions.tsv", {"H3K27ac": "ag_h3k27ac_mean_score"}),
        "RNA_deletions": ("FigureS5D_RNA_deletions.tsv", {
            "three_gene_mean": "ag_rna_3gene_mean_percent_change", "SORT1": "ag_rna_percent_change_SORT1",
            "PSRC1": "ag_rna_percent_change_PSRC1", "CELSR2": "ag_rna_percent_change_CELSR2"}),
    }
    for analysis, (filename, endpoints) in s5_specs.items():
        path = f"outputs/source_data/{filename}"
        data = pd.read_csv(ROOT / path, sep="\t")
        for regime, part in data.groupby("model"):
            for endpoint, column in endpoints.items():
                valid = part["kircher_mean_log2_effect"].notna() & part[column].notna()
                add(f"figs5_{analysis.lower()}_{regime.lower()}_{endpoint.lower()}_pearson", "figs5_kircher_model_exposure",
                    f"Kircher versus AlphaGenome Pearson r, {analysis}, {endpoint}, {regime}",
                    pearsonr(part.loc[valid, "kircher_mean_log2_effect"], part.loc[valid, column]).statistic,
                    "Pearson r", regime, f"n={int(valid.sum())} edits", path)
                add(f"figs5_{analysis.lower()}_{regime.lower()}_{endpoint.lower()}_spearman", "figs5_kircher_model_exposure",
                    f"Kircher versus AlphaGenome Spearman rho, {analysis}, {endpoint}, {regime}",
                    spearmanr(part.loc[valid, "kircher_mean_log2_effect"], part.loc[valid, column]).statistic,
                    "Spearman rho", regime, f"n={int(valid.sum())} edits", path)

    # Figure S6: per-element matched held-out accessibility benchmarks.
    s6_path = "outputs/source_data/FigureS6_kircher_other_loci/FigureS6C_accessibility_correlation_model_comparison_source.tsv"
    s6 = pd.read_csv(ROOT / s6_path, sep="\t")
    for _, row in s6.iterrows():
        key = f"{row['element']}_{row['modality']}_{row['model']}".lower().replace(" ", "_")
        add(f"figs6_{key}_rho", "figs6_kircher_matched_heldout",
            f"Element-level Spearman rho, {row['element']}, {row['modality']}, {row['model']}", row["spearman_rho"],
            "Spearman rho", row["model"], f"n={int(row['n'])} substitutions", s6_path,
            f"Assigned held-out fold: {row['model_version']}" if pd.notna(row["model_version"]) else "Default ALL_FOLDS model.")

    # Supplementary quantitative controls used explicitly in the text.
    s9a_path = "outputs/source_data/FigureS9_module_transfer_controls/S9A_values.csv"
    s9a = pd.read_csv(ROOT / s9a_path)
    add("figs9a_pearson", "figs9_module_transfer_controls", "HPA versus AlphaGenome native RNA Pearson r",
        pearsonr(s9a["hpa_log10"], s9a["ag_log10"]).statistic, "Pearson r", "ALL_FOLDS",
        f"n={len(s9a)} genes", s9a_path)
    add("figs9a_spearman", "figs9_module_transfer_controls", "HPA versus AlphaGenome native RNA Spearman rho",
        spearmanr(s9a["hpa_log10"], s9a["ag_log10"]).statistic, "Spearman rho", "ALL_FOLDS",
        f"n={len(s9a)} genes", s9a_path)
    for panel, filename, ycol in (("D", "S9D_values.csv", "T_minus_G"),
                                  ("E", "S9E_values.csv", "log2_T_over_G")):
        path = f"outputs/source_data/FigureS9_module_transfer_controls/{filename}"
        data = pd.read_csv(ROOT / path)
        add(f"figs9{panel.lower()}_spearman", "figs9_module_transfer_controls",
            f"Native RNA versus {ycol} Spearman rho", spearmanr(data["native"], data[ycol]).statistic,
            "Spearman rho", "ALL_FOLDS", f"n={len(data)} recipients", path)

    s10_path = "outputs/source_data/FigureS10_distal_and_tissue_controls/fold0_tissue_matrix.tsv"
    fold0 = pd.read_csv(ROOT / s10_path, sep="\t").set_index("context")
    aligned_all = tissue.loc[fold0.index, fold0.columns].to_numpy().ravel()
    aligned_f0 = fold0.to_numpy().ravel()
    add("figs10b_spearman", "figs10_tissue_replication", "ALL_FOLDS versus FOLD_0 tissue-gene Spearman rho",
        spearmanr(aligned_all, aligned_f0).statistic, "Spearman rho", "ALL_FOLDS;FOLD_0",
        f"n={aligned_all.size} tissue-gene cells", s10_path)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} verified headline results to {OUT}")


if __name__ == "__main__":
    main()
