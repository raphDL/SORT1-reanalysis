# Manifest audit notes

`MANIFEST.tsv` is the authoritative panel crosswalk for the release. It is
assembled from the latest accepted panel decisions, not from legacy folder
numbers alone.

## Status meanings

- `final` — current scientific content and panel assignment are accepted.
- `author-layout` — the panel is an author-composed schematic or adapted
  graphical abstract; the editable Inkscape source must be added manually.
- `exclude` — retained in the working archive but not selected for release.

## Known numbering hazard

Several working folders reuse `S5` and `S6` for different generations. The
release therefore identifies panels by scientific content and canonical asset,
not by the legacy directory name. The final numbering was reconciled against
the author-approved manuscript package frozen on 2026-08-01.

The inferred current supplementary order is:

1. contact controls and cross-context contact predictions;
2. eQTL/tagging and FOLD_0 variant-prioritization diagnostics;
3. liver caQTL locus analysis;
4. human CRISPR repair weighting and repair-junction reconstruction;
5. SORT1 Kircher model comparison;
6. other-locus Kircher model comparison;
7. native-sequence substitution and motif audit;
8. scramble-boundary and component-necessity robustness;
9. module-transfer controls;
10. distal-contact and tissue-context controls.

## Panels assembled outside Python

Adapted graphical abstracts and simple design cartoons were composed in
Inkscape. They must be deposited as editable SVGs and mapped here before the
release is declared complete. A raster image embedded in the manuscript is not
sufficient as the archival source.

## Release rule

No file is copied from the 38-GB working archive merely because it exists.
Inclusion requires one of the following:

1. it is a published panel or its compact source table;
2. it is needed to regenerate a reported number;
3. it records a frozen decision, model assignment, provenance check, or
   prediction checksum;
4. it is a small curated input that cannot be fetched reproducibly.

## Data-copying pass (2026-08-01)

Compact source tables for every non-`author-layout`, non-`exclude` panel have
been copied into `outputs/source_data/` with SHA-256 checksums recorded
alongside the copy log. Tables that are not truly compact (>15 MB; not raw
prediction tensors, but full per-position/per-variant scans) were **not**
copied into git. Their SHA-256 digests, sizes, and legacy paths are recorded
in `outputs/run_manifests/zenodo_pending_large_outputs.tsv` instead, pending
an actual Zenodo deposit and DOI:

- Fig3A: `figure2A_full_region_snv_position_summary.tsv` (171 MB) and
  `figure2A_full_region_snv_synergy_liver_best_alt_by_position.tsv` (18.5 MB)
- Fig4H: `liver_cd14_monocyte_tcell_full_region_synergy_same_axis_130mm.tsv` (17.3 MB)
- Fig4J: `dense_100kb_8tf_3tracks_variant_summary.tsv` (528 MB)

Fig1B originally had **no compact source table anywhere in the working
archive**: its legacy script (`panel_tracks/panel_a_sort1_locus_rnaseq_cebp.py`)
plotted AlphaGenome `TrackData` directly from a live API call. A one-time
authorized `predict_variant` regeneration (524,288 bp, rs12740374,
ALL_FOLDS) has now populated `outputs/source_data/Figure1B_locus_tracks/`.
The release exporter reads its credential only from the environment and does
not carry forward the hardcoded API-key fallback found in the legacy script.

## Script-porting status (2026-08-01)

`figures/` is being populated per figure. Release renderers read only frozen
files under `outputs/source_data/`; they never call AlphaGenome or download
public data. The frozen author-approved assets remain authoritative for final
multi-panel composition.

**Ported and test-rendered:** Fig1C–F (`figures/fig1.py`), Fig2B/C/E/F
(`figures/fig2.py`), Fig3E–G (`figures/fig3.py`) and FigS5B/C
(`figures/figS5.py`). These were inspected for scientific and visual
equivalence. Exact byte identity with the legacy export is not claimed unless
the manifest explicitly says so, because plotting-library metadata and object
IDs differ.

**Not yet ported** (compact source tables are already copied except where
explicitly stated):

| Panel | Legacy script | Why not ported this pass |
|---|---|---|
| S1A/S1C/S1D | `panel_contact/run_observed_hic_validation.py` | Supplementary panels remain represented by their frozen canonical SVGs. |
| S1B | `supplementary_figures/build_figure1_supplementary_panels.py` | Shared multi-panel supplementary-figure builder. |
| S2A-C | same file as S1B | Same builder, different function. |
| 3A | `figure3_restructured/panel_A_regional_coordinated_ism/make_panel_A.py` | Depends on a shared helper module (`make_compact_global_liver_ism_panel.py`) plus the two Zenodo-pending large tables above and a GENCODE feather cache; not self-contained with committed release data. |
| 3B | `figure3_restructured/panel_B_base_substitution_501bp/run_native_locus_501bp_ism.py` | 860 lines, combined analysis+plot. |
| 3C | `figure3_restructured/panel_C_motif_family_disruption/make_panel_C.py` | Requires a live genome FASTA (`pysam`) + JASPAR PFM scan, not a pure replot; belongs in `analysis/`, not `figures/`. |

Figure 4 and supplementary panels not listed above retain canonical SVGs,
compact source tables, hashes and exact legacy-script pointers. Further
renderer extraction is useful but is not a numerical-provenance blocker.
