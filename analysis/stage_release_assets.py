#!/usr/bin/env python3
"""Stage unambiguous panel assets from the frozen working archive.

This migration helper copies only rows representing a single panel and having
an explicit, existing legacy asset.  Composite supplementary rows (for
example ``S2A-C``) and author-layout placeholders are deliberately skipped:
their final assembly must be reconciled against the submission PDF.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PORTED_RENDERERS = {"3E", "3F", "3G", "S5B", "S5C"}
COMPOSITE_COMPONENTS = {
    "S2A-C": {
        "FigureS2A.svg": "investigation/SORT1_cholesterol_musunuru_2010/report/supplementary_figures/panel_figure1_fold0_variant_prioritization/figureS1A_gtex_vs_rs127_ld_tagging_all_folds.svg",
        "FigureS2B.svg": "investigation/SORT1_cholesterol_musunuru_2010/report/supplementary_figures/panel_figure1_fold0_variant_prioritization/figureS1B_gtex_vs_rs127_ld_tagging_fold_0.svg",
        "FigureS2C.svg": "investigation/SORT1_cholesterol_musunuru_2010/report/supplementary_figures/panel_figure1_fold0_variant_prioritization/figureS1C_all_folds_vs_fold0.svg",
    },
    "S3A-G": {
        f"FigureS3{letter}.svg": f"investigation/SORT1_cholesterol_musunuru_2010/report/supplementary_figures/panel_s2_caqtl_locus/figureS2{letter}_{stem}.svg"
        for letter, stem in (
            ("A", "observed_caqtl"), ("B", "all_folds_direct_atac"),
            ("C", "all_folds_tagging"), ("D", "observed_vs_all_folds_tagging"),
            ("E", "fold0_direct_atac"), ("F", "fold0_tagging"),
            ("G", "observed_vs_fold0_tagging"),
        )
    },
    "S6A-C": {
        "FigureS6A.svg": "investigation/SORT1_cholesterol_musunuru_2010/report/supplementary_figures/panel_s6_kircher_other_loci_model_comparison/FigureS6A_DNase_ALL_FOLDS_vs_matched_heldout.svg",
        "FigureS6B.svg": "investigation/SORT1_cholesterol_musunuru_2010/report/supplementary_figures/panel_s6_kircher_other_loci_model_comparison/FigureS6B_ATAC_ALL_FOLDS_vs_matched_heldout.svg",
        "FigureS6C.svg": "investigation/SORT1_cholesterol_musunuru_2010/report/supplementary_figures/panel_s6_kircher_other_loci_model_comparison/FigureS6C_accessibility_correlation_model_comparison.svg",
    },
    "S8A-E": {
        "FigureS8A.svg": "investigation/SORT1_cholesterol_musunuru_2010/report/supplementary_figures/panel_s9_scramble_boundary_robustness/FigureS9A_three_gene_directional_recovery.svg",
        "FigureS8B.svg": "investigation/SORT1_cholesterol_musunuru_2010/report/supplementary_figures/panel_s9_scramble_boundary_robustness/FigureS9B_three_gene_arm_necessity_ALL_FOLDS_FOLD0.svg",
        "FigureS8C.svg": "investigation/SORT1_cholesterol_musunuru_2010/report/supplementary_figures/panel_s9_scramble_boundary_robustness/FigureS9C_complete_2816_window_heatmap.svg",
        "FigureS8D.svg": "investigation/SORT1_cholesterol_musunuru_2010/report/supplementary_figures/panel_s9_scramble_boundary_robustness/FigureS9D_seed_level_313_314_315bp_thresholds.svg",
        "FigureS8E.svg": "investigation/SORT1_cholesterol_musunuru_2010/report/supplementary_figures/panel_s9_scramble_boundary_robustness/FigureS9E_final_315bp_three_gene_seed_distributions.svg",
    },
    "S9A-E": {
        "FigureS9A.svg": "investigation/SORT1_cholesterol_musunuru_2010/report/supplementary_figures/panel_s10_module_transfer_controls/panel_A_HPA_vs_AlphaGenome_native_RNA/FigureS10A_HPA_liver_vs_AlphaGenome_native_RNA.svg",
        "FigureS9B.svg": "investigation/SORT1_cholesterol_musunuru_2010/report/supplementary_figures/panel_s10_module_transfer_controls/panel_B_distance_distributions/FigureS10B_T_minus_G_distribution_by_distance.svg",
        "FigureS9C.svg": "investigation/SORT1_cholesterol_musunuru_2010/report/supplementary_figures/panel_s10_module_transfer_controls/panel_C_paired_T_minus_G/FigureS10C_paired_T_minus_G_by_HPA_cohort.svg",
        "FigureS9D.svg": "investigation/SORT1_cholesterol_musunuru_2010/report/supplementary_figures/panel_s10_module_transfer_controls/panel_D_native_vs_absolute_effect/FigureS10D_native_RNA_vs_absolute_T_minus_G.svg",
        "FigureS9E.svg": "investigation/SORT1_cholesterol_musunuru_2010/report/supplementary_figures/panel_s10_module_transfer_controls/panel_E_native_vs_relative_effect/FigureS10E_native_RNA_vs_relative_T_over_G.svg",
    },
    "S10A-B": {
        "FigureS10A.svg": "investigation/SORT1_cholesterol_musunuru_2010/report/supplementary_figures/panel_s11_distal_and_context_controls/FigureS11A_contact_associated_RNA_by_50kb_distance.svg",
        "FigureS10B.svg": "investigation/SORT1_cholesterol_musunuru_2010/report/supplementary_figures/panel_s11_distal_and_context_controls/FigureS11B_rs12740374_tissue_replication_ALL_FOLDS_FOLD0.svg",
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def is_single_panel(panel: str) -> bool:
    return "-" not in panel and len(panel) >= 2 and panel[-1].isalpha()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace-root", type=Path, required=True)
    args = parser.parse_args()

    rows = list(csv.DictReader((ROOT / "MANIFEST.tsv").open(), delimiter="\t"))
    records: list[dict[str, str]] = []
    for row in rows:
        panel = row["panel"]
        legacy = row["legacy_asset"]
        release = row["release_asset"]
        if row["status"] == "exclude" or panel in PORTED_RENDERERS or not is_single_panel(panel):
            continue
        if not legacy or legacy.startswith(("TO_", "PENDING_")) or not release:
            continue
        source = args.workspace_root / legacy
        destination = ROOT / release
        if not source.is_file():
            continue
        if source.suffix.lower() != destination.suffix.lower():
            # Never place PDF bytes in a file named .svg (or vice versa).
            # The manifest must name the release asset's real format.
            continue
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
        records.append(
            {
                "panel": panel,
                "release_asset": release,
                "sha256": sha256(destination),
                "legacy_asset": legacy,
                "operation": "byte-identical copy",
            }
        )

    for panel, components in COMPOSITE_COMPONENTS.items():
        destination_dir = ROOT / f"figures/rendered/Figure{panel.split('A-')[0]}"
        destination_dir.mkdir(parents=True, exist_ok=True)
        for filename, legacy in components.items():
            source = args.workspace_root / legacy
            if not source.is_file():
                raise FileNotFoundError(source)
            destination = destination_dir / filename
            shutil.copyfile(source, destination)
            records.append({
                "panel": panel,
                "release_asset": str(destination.relative_to(ROOT)),
                "sha256": sha256(destination),
                "legacy_asset": legacy,
                "operation": "byte-identical component copy",
            })

    # S7 is a single accepted audit panel whose working filename still says S8.
    s7_source_rel = "investigation/SORT1_cholesterol_musunuru_2010/report/supplementary_figures/panel_s8_native_sequence_sensitivity/FigureS8B_ALL_FOLDS_vs_FOLD0_hotspot_audit.svg"
    s7_source = args.workspace_root / s7_source_rel
    s7_destination = ROOT / "figures/rendered/FigureS7.svg"
    shutil.copyfile(s7_source, s7_destination)
    records.append({
        "panel": "S7", "release_asset": str(s7_destination.relative_to(ROOT)),
        "sha256": sha256(s7_destination), "legacy_asset": s7_source_rel,
        "operation": "byte-identical copy with final-number rename",
    })

    out = ROOT / "outputs/run_manifests/release_asset_provenance.tsv"
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("panel", "release_asset", "sha256", "legacy_asset", "operation"),
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(records)
    print(f"Staged {len(records)} unambiguous panel assets; wrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
