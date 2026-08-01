#!/usr/bin/env python3
"""Add explicit output, track, scorer, window, and effect metadata per analysis."""

from __future__ import annotations

import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / "outputs" / "run_manifests" / "panel_analysis_parameters.tsv"


def main() -> None:
    with PATH.open(encoding="utf-8", newline="") as handle:
        base = list(csv.DictReader(handle, delimiter="\t"))

    # Values are intentionally analysis-level rather than inferred from figure labels.
    # Semicolon-separated entries describe multiple outputs used by one panel.
    meta: dict[str, tuple[str, str, str, str, str]] = {
        "1B": ("RNA;CHIP_TF;DNASE;ATAC", "HepG2 CEBPA; liver UBERON:0002107",
               "Displayed tracks selected by assay and ontology metadata",
               "Raw sequence prediction (no variant scorer)", "Displayed 70-kb locus interval"),
        "1C": ("RNA_SEQ", "liver UBERON:0002107",
               "Gene-strand-compatible indices: CELSR2 233/234/627; PSRC1 and SORT1 504/505/627",
               "GeneMaskLFCScorer", "Exon-masked gene score; mean of 3 compatible tracks per gene"),
        "1E": ("CONTACT_MAPS", "HepG2 EFO:0001187", "All HepG2 contact-map tracks; averaged",
               "Raw FOLD_0 contact-map output", "2,048-bp bins; 3x3 anchor/promoter mean"),
        "1F": ("CONTACT_MAPS", "HepG2 EFO:0001187; experimental 4DNFICSTCJQZ",
               "Predicted HepG2 tracks averaged; observed KR-balanced O/E map",
               "Same-distance percentile", "2/2.048-kb maps; 3x3 anchor/promoter mean"),
        "2B": ("RNA_SEQ", "liver UBERON:0002107",
               "Primary target index 4: pooled child/adult liver poly(A)+ unstranded RNA-seq",
               "Raw sequence prediction (no variant scorer)", "Mean signal within +/-2 kb of each TSS"),
        "2C": ("RNA_SEQ", "liver UBERON:0002107",
               "Primary target index 4: pooled child/adult liver poly(A)+ unstranded RNA-seq",
               "Raw sequence prediction (no variant scorer)", "Mean signal within +/-2 kb of each TSS"),
        "2E": ("RNA_SEQ;ATAC;CHIP_HISTONE", "HepG2 EFO:0001187",
               "Gene-strand-compatible RNA; matched HepG2 ATAC; H3K27ac-selected histone tracks",
               "Recommended RNA/ATAC/histone variant scorers", "RNA gene score; mean signed ATAC/H3K27ac score"),
        "2F": ("DNASE;ATAC", "Experimental cell-line matched",
               "DNase for all 6 elements; matched ATAC for F9/LDLR/PKLR/SORT1",
               "Center-mask DIFF_LOG2_SUM", "16,384-bp native-locus substitution interval"),
        "3A": ("RNA_SEQ", "liver UBERON:0002107", "All liver tracks; averaged",
               "Raw sequence prediction (no variant scorer)", "Mean signal within +/-2 kb of 3 TSSs"),
        "3B": ("RNA_SEQ", "liver UBERON:0002107", "All liver tracks; averaged",
               "Raw sequence prediction (no variant scorer)", "Mean signal within +/-2 kb of each TSS; equal 3-gene mean"),
        "3C": ("derived PWM", "JASPAR 2024 CORE vertebrates", "Frozen H1-H6 and C/EBP windows",
               "PWM compatibility ranking", "Derived from 3B maximum-loss substitutions"),
        "3E": ("RNA_SEQ", "liver UBERON:0002107", "All liver tracks; averaged",
               "Raw sequence prediction (no variant scorer)", "TSS +/-2 kb; allele-effect retention; 8 seeds"),
        "3F": ("RNA_SEQ", "liver UBERON:0002107", "All liver tracks; averaged",
               "Raw sequence prediction (no variant scorer)", "SORT1 TSS +/-2 kb; allele-effect retention; 8 seeds"),
        "3G": ("RNA_SEQ", "liver UBERON:0002107", "All liver tracks; averaged",
               "Raw sequence prediction (no variant scorer)", "TSS +/-2 kb; equal 3-gene retention mean; 8 seeds"),
        "4B": ("RNA_SEQ", "liver UBERON:0002107", "Gene-strand-compatible liver tracks; averaged",
               "Raw sequence prediction (no variant scorer)", "Mean signal within +/-2 kb of recipient TSS"),
        "4C": ("RNA_SEQ", "liver UBERON:0002107", "Gene-strand-compatible liver tracks; averaged",
               "Raw sequence prediction (no variant scorer)", "Mean signal within +/-2 kb of recipient TSS"),
        "4E": ("RNA_SEQ", "HepG2 EFO:0001187", "Promoter-compatible HepG2 RNA output",
               "Raw sequence prediction (no variant scorer)", "4,001-bp recipient-TSS window; high-minus-low T-minus-G interaction"),
        "4F": ("RNA_SEQ", "HepG2 EFO:0001187", "Promoter-compatible HepG2 RNA output",
               "Raw sequence prediction (no variant scorer)", "Same endpoint as 4E; within-distance contact-contrast quintiles"),
        "4G": ("RNA_SEQ", "7 tissue ontologies", "Gene-strand-compatible tracks averaged within ontology",
               "Recommended RNA variant scorer", "Gene-level G-to-T effect"),
        "4H": ("RNA_SEQ", "liver; CD8+ memory T cell; CD14+ monocyte",
               "Context-matched RNA tracks", "Raw sequence prediction (no variant scorer)",
               "TSS +/-2 kb; minimum positive effect shared by all 3 genes"),
        "4J": ("RNA_SEQ", "liver; CD8+ memory T cell; CD14+ monocyte",
               "Context-matched RNA tracks", "Raw sequence prediction (no variant scorer)",
               "TSS +/-2 kb; signed 3-gene sum; best signed orientation"),
        "S1B": ("CONTACT_MAPS", "HepG2 EFO:0001187; experimental 4DNFICSTCJQZ",
                 "Same tracks/maps as 1E-1F", "Same-distance percentile", "Single TSS-containing bin"),
        "S1C": ("CONTACT_MAPS", "HepG2 EFO:0001187", "HepG2 contact tracks",
                 "ContactMapScorer plus promoter-specific O/E delta", "2,048-bp promoter bins; local 100-SNV control"),
        "S1D": ("CONTACT_MAPS", "12 available model ontologies", "All 28 tracks, averaged within ontology",
                 "Raw reference contact-map output", "2,048-bp bins; 3-bin promoter summary; same-distance percentile"),
        "S2A-C": ("RNA_SEQ", "liver UBERON:0002107",
                   "Same gene-strand-compatible indices as 1C", "GeneMaskLFCScorer",
                   "Exon-masked gene score; 111 direct variants; no LD term in panel C"),
        "S3A-G": ("ATAC", "liver UBERON:0002107; HepG2 EFO:0001187",
                   "All liver and HepG2 ATAC outputs; summed", "Center-mask DIFF_LOG2_SUM",
                   "Variant-level direct or rs12740374-tagging score"),
        "S4A": ("experimental read counts", "primary human hepatocytes", "Wang et al. human reads",
                 "No AlphaGenome scoring in panel A", "Human repair-depth and product coverage summary"),
        "S4B": ("RNA_SEQ", "liver UBERON:0002107",
                 "Target index 4; same cached predictions as human repair analysis", "Raw sequence prediction",
                 "TSS +/-2 kb; product effect multiplied by conditional human indel frequency"),
        "S4C": ("RNA_SEQ", "liver UBERON:0002107", "Target index 4", "Raw sequence prediction",
                 "TSS +/-2 kb; selected deletion geometries"),
        "S4D": ("RNA_SEQ", "liver UBERON:0002107", "Target index 4", "Raw sequence prediction",
                 "TSS +/-2 kb; selected deletion geometries grouped by junction status"),
        "S5A": ("RNA_SEQ", "HepG2 EFO:0001187", "Gene-strand-compatible RNA tracks",
                 "Recommended RNA variant scorer", "Construct-averaged MPRA substitutions; 3 genes and mean"),
        "S5B": ("ATAC;CHIP_HISTONE", "HepG2 EFO:0001187", "Matched ATAC and H3K27ac tracks",
                 "Recommended ATAC/histone variant scorers", "Construct-averaged MPRA substitutions"),
        "S5C": ("RNA_SEQ", "HepG2 EFO:0001187", "Gene-strand-compatible RNA tracks",
                 "Recommended RNA scorer on reconstructed one-base deletions", "Construct-averaged MPRA deletions"),
        "S5D": ("RNA_SEQ", "HepG2 EFO:0001187", "Gene-strand-compatible RNA tracks",
                 "Recommended RNA scorer on reconstructed one-base deletions", "Legacy grouped S5 range; final figure has A-C"),
        "S6A-C": ("DNASE;ATAC", "Experimental cell-line matched", "Same per-element tracks as 2F",
                   "Center-mask DIFF_LOG2_SUM", "Identical edits under ALL_FOLDS and matched held-out fold"),
        "S7": ("RNA_SEQ", "liver UBERON:0002107", "All liver tracks; averaged",
               "Raw sequence prediction", "Combined-window loss at 7 frozen windows; 3 TSSs +/-2 kb"),
        "S8A-E": ("RNA_SEQ", "liver UBERON:0002107", "All liver tracks; averaged",
                   "Raw sequence prediction", "TSS +/-2 kb; scramble retention; 8 seeds"),
        "S9A-E": ("RNA_SEQ", "liver UBERON:0002107", "Gene-strand-compatible liver tracks; averaged",
                   "Raw sequence prediction", "Recipient TSS +/-2 kb; native/transfer controls"),
        "S10A-B": ("RNA_SEQ", "HepG2; 7 tissue ontologies",
                    "S10A HepG2 RNA; S10B gene-strand-compatible tissue tracks",
                    "Raw sequence prediction; recommended RNA variant scorer",
                    "S10A contact-associated T-minus-G; S10B gene-level G-to-T effects"),
    }

    out = []
    for row in base:
        panel = row["panel"]
        if panel not in meta:
            raise KeyError(f"No explicit analysis metadata for {panel}")
        modality, context, tracks, scorer, aggregation = meta[panel]
        out.append({
            "panel": panel,
            "model_regime": row["model_regime"],
            "prediction_context_bp": row["prediction_context_bp"],
            "output_modality": modality,
            "ontology_or_biosample": context,
            "track_selection": tracks,
            "scorer_or_measure": scorer,
            "aggregation_window_or_endpoint": aggregation,
            "output_summary": row["output_summary"],
            "parameter_source": row["parameter_source"],
            "fasta_source": row["fasta_source"].replace(
                "NEEDS_MANUAL_CHECK — rendering script itself performs no sequence fetch, but plots output of the local_fasta-based Wang analysis (2B); confirm no independent construction step",
                "derived_from_local_fasta_predictions — renderer uses cached outputs from the same Wang reconstruction pipeline as 2B",
            ),
        })

    with PATH.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(out)
    print(f"Wrote explicit analysis metadata for {len(out)} panel groups to {PATH}")


if __name__ == "__main__":
    main()
