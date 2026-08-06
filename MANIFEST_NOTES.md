# Manifest notes

`MANIFEST.tsv` is the panel crosswalk for the release: for every figure
panel, the source table, render script, model regime, and working-archive
origin.

## Status values

- `final` — scientific content and panel assignment are settled.
- `author-layout` — an author-composed schematic or adapted graphical
  abstract; the archival source is an editable Inkscape SVG, not a script.
- `exclude` — present in the working archive but not part of the release.

## Panel numbering

Panels are identified by scientific content and canonical asset rather
than legacy directory name, since some working-archive folder names were
reused across analysis generations. Supplementary panels are numbered
S1–S10, matching the manuscript.

## Release inclusion rule

A file from the working archive is included only if it is:

1. a published panel or its compact source table;
2. needed to regenerate a reported number;
3. a frozen decision, model assignment, provenance check, or prediction
   checksum; or
4. a small curated input that cannot be fetched reproducibly.

## Large derived tables (deposited on Zenodo)

Compact source tables (<15 MB) for every non-`author-layout`,
non-`exclude` panel are in `outputs/source_data/` with SHA-256 checksums.
Full per-position/per-variant/per-track scan tables too large to commit
are archived at Zenodo: [10.5281/zenodo.21820090](https://doi.org/10.5281/zenodo.21820090).
Checksums, sizes, and legacy paths for each file are in
`outputs/run_manifests/zenodo_pending_large_outputs.tsv`:

- Fig3A: `figure2A_full_region_snv_position_summary.tsv` (171 MB) and
  `figure2A_full_region_snv_synergy_liver_best_alt_by_position.tsv` (18.5 MB)
- Fig3F: 1bp boundary-window scan, all folds (3.3 MB zipped)
- Fig4E/4F: full per-track/per-state RNA-seq/ATAC/H3K27ac predictions for
  the 315bp T/G transfer benchmark (`track_level.tsv` 31 MB,
  `predictions_raw.tsv` 9 MB, `activity_raw.tsv` 6.3 MB)
- Fig4H: `liver_cd14_monocyte_tcell_full_region_synergy_same_axis_130mm.tsv` (17.3 MB)
- Fig4J: `dense_100kb_8tf_3tracks_variant_summary.tsv` (528 MB)
- FigS8D: 2D scramble-boundary coordinate-grid and minimal-core sweeps,
  all folds (10.9 MB + 4.6 MB zipped)

Figure 1B's legacy script plotted AlphaGenome track data directly from a
live API call with no cached source table; `outputs/source_data/
Figure1B_locus_tracks/` was regenerated via `predict_variant` (524,288 bp,
rs12740374, `ALL_FOLDS`).

## Deterministic figure-rendering scripts

`figures/fig*.py` render published panels from `outputs/source_data/`
only — no AlphaGenome calls, no downloads. Ported so far: Fig1C–F
(`fig1.py`), Fig2B/C/E/F (`fig2.py`), Fig3E–G (`fig3.py`), FigS5B/C
(`figS5.py`). Not yet ported (their compact source tables are already in
the release; only the plotting script is pending):

| Panel | Legacy script | Note |
|---|---|---|
| S1A/S1C/S1D | `panel_contact/run_observed_hic_validation.py` | Represented by frozen canonical SVGs for now. |
| S1B, S2A–C | `build_figure1_supplementary_panels.py` | Shared multi-panel builder. |
| 3A | `panel_A_regional_coordinated_ism/make_panel_A.py` | Depends on the two Zenodo-pending Fig3A tables above plus a GENCODE feather cache. |
| 3B | `panel_B_base_substitution_501bp/run_native_locus_501bp_ism.py` | Combined analysis+plot script. |
| 3C | `panel_C_motif_family_disruption/make_panel_C.py` | Needs a live genome FASTA + JASPAR scan, not a pure replot. Ported to `reproduction/figure3.py::run_fig3c` as part of the clean-room analysis runner instead. |

Panels not listed above (including all of Figure 4's computational
panels) retain canonical SVGs, compact source tables, checksums, and
legacy-script pointers in `MANIFEST.tsv`.
