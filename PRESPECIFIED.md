# Frozen and prespecified analyses

This file is the root-level index for decisions that were frozen before the
corresponding validation or sensitivity result was inspected. The exact
original artifacts, hashes, and dates preserved for this release are in
`outputs/run_manifests/`.

## Artifacts to preserve

| Analysis | Original artifact |
|---|---|
| Multi-element Kircher benchmark | `panel_kircher_multielement/PRESPECIFIED_PLAN.md` |
| Matched held-out model assignment | `panel_kircher_multielement/heldout_matched_folds/FROZEN_MODEL_ASSIGNMENT.json` |
| Figure 4 motif candidates and matched shuffles | `panel_TF_insertion/PRESPECIFIED_VALIDATION_PLAN.md`; `prespecified_*` TSV/JSON files |
| Genome-wide caQTL benchmark | `panel_currin_genomewide_benchmark/PRESPECIFIED_PLAN.md`; frozen manifest and SHA-256 |
| FOLD_0 tissue replication | `panel_fold0_tissue_replication/PRESPECIFIED_PLAN.md` |
| FOLD_0 native hotspot audit | `panel_fold0_native_hotspots/PRESPECIFIED_PLAN.md`; frozen candidates and SHA-256 |
| Distal Hi-C transfer track selection | `panel_distal_hic_transfer/FROZEN_TRACK_SELECTION_V2.json` |
| Distal Hi-C robustness controls | `panel_distal_hic_transfer/results/chr1_distal_315_transfer_robustness/frozen_control_metadata.json` |
| Frozen Figure 4E-F statistics | `panel_distal_hic_transfer/main_figure4_panels/frozen_v1/` |

## Release requirements

For each frozen analysis, the public copy must preserve:

- the original timestamp and content;
- the candidate or sample manifest;
- the relevant SHA-256 digest;
- the model regime and AlphaGenome code revision;
- whether the criterion passed, failed, or was underpowered;
- any post-freeze correction, with the original retained and the correction
  described rather than overwritten.

The existence of this index is not evidence that every analysis was genuinely
prespecified. Only analyses with contemporaneous source artifacts are labeled
as such in the manuscript.
