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
| S1A–D | `panel_contact/run_observed_hic_validation.py`, `build_figure1_supplementary_panels.py`, `panel_contact/run_contact_architecture_analysis.py` | Ported to `reproduction/figureS1.py::run_figs1a`–`run_figs1d` as part of the clean-room analysis runner instead. |
| S2A–C | `build_figure1_supplementary_panels.py` | Ported to `reproduction/figureS2.py::run_figs2a`–`run_figs2c` as part of the clean-room analysis runner instead. See "Known temporal drift" below for a caveat specific to S2C. |
| S3A–G | `panel_caQTL/run_currin_score_249_ag.py`, `panel_caQTL/run_currin_caqtl_atac_analysis.py`, `build_figure_s2.py` | Ported to `reproduction/figureS3.py::run_figs3` as part of the clean-room analysis runner instead. Needs `--currin-variants`/`--currin-peakset` (see "External data requiring manual supply" below). |
| S4A–D | `panel_s4_wang_repair_weighting/make_figure_s4_weighted_panels.py`, `panel_wang_et_al_indels/plot_human_frequency_weighted_heatmap.py`, `archive_s4_crispr_motif_and_insertions_before_split/make_figure_s4_panels.py`, `synthetic_cutsite_scan/run_deletion_xy_full_grid.py` | Ported to `reproduction/figureS4.py::run_figs4a`–`run_figs4d`. Makes no AlphaGenome calls of its own -- S4B/C/D reuse Figure 2B/2C's already-scored sequences (verified: the 24 human indels and the 30 selected deletion designs are both within Figure 2B/2C's own scan). S4A's total/unedited read counts are frozen (no raw-sequencing source in the redistributed Wang spreadsheet); see the module docstring. |
| 3A | `panel_A_regional_coordinated_ism/make_panel_A.py` | Depends on the two Zenodo-pending Fig3A tables above plus a GENCODE feather cache. |
| 3B | `panel_B_base_substitution_501bp/run_native_locus_501bp_ism.py` | Combined analysis+plot script. |
| 3C | `panel_C_motif_family_disruption/make_panel_C.py` | Needs a live genome FASTA + JASPAR scan, not a pure replot. Ported to `reproduction/figure3.py::run_fig3c` as part of the clean-room analysis runner instead. |

Panels not listed above (including all of Figure 4's computational
panels) retain canonical SVGs, compact source tables, checksums, and
legacy-script pointers in `MANIFEST.tsv`.

## External data requiring manual supply

Two inputs have no traceable derivation script from their public release to
the small local excerpt actually used, so `reproduce.py` cannot fetch them
automatically; they must be supplied via a CLI flag and are only checksum-
validated, never committed to this repository (`data/README.md`: "No public
dataset should be committed merely for convenience"):

- **Fig2B** (`--wang-xls`): the original Wang et al. 2018 supplementary
  spreadsheet (`atv310103_ds.xls`); the publisher download is bot-gated.
- **FigS3** (`--currin-variants`, `--currin-peakset`): a SORT1-locus excerpt
  of the Currin et al. 2025 liver caQTL summary statistics and its 28-peak
  "coordinated set 539" definition (Zenodo record 15025748 / GEO GSE264684;
  accession and checksums in `data/SOURCES.tsv`). The full Zenodo release is
  genome-wide; reducing it to this locus and peak set was done once, by
  hand, outside this codebase.

## Known temporal drift in the AlphaGenome backend (observed 2026-08-06)

The AlphaGenome API is a live service; its deployed model is not
guaranteed to return byte-identical predictions across calendar time, even
for a fixed model version (`ALL_FOLDS`/`FOLD_0`) and a fixed input. This
was directly observed while reproducing Figure 1C/S2's 111-variant
liver exon-mask RNA scan:

- A clean-room run from 2026-08-02 matched the frozen `outputs/source_data/`
  reference to ~1e-12 (floating-point noise) for every one of the 999
  variant x gene x track predictions underlying Figure 1C.
- An otherwise-identical clean-room run from 2026-08-06 (4 days later)
  diverged from that same reference by up to 0.044 at the individual-track
  level and up to ~0.006 in the gene-level aggregate (Pearson r 0.98–0.999).
  The same 111 liver tracks were returned in both runs, so this is
  per-prediction numerical drift, not a track-set change.
- Two independent fresh FOLD_0 scans taken ~2h apart on 2026-08-06 were
  bit-identical (max abs diff 0.0 across all 111 variants x 3 genes),
  confirming the 2026-08-06 predictions are themselves stable/deterministic
  -- the divergence is a real shift in the deployed backend between
  2026-08-02 and 2026-08-06, not request-to-request sampling noise.

**Consequence for S2C specifically:** the manuscript states "rs12740374
remained the variant with the largest absolute predicted RNA effect for
SORT1, CELSR2, and PSRC1" under held-out FOLD_0. This holds for CELSR2 and
PSRC1 in both the 2026-08-02 and 2026-08-06 runs. For SORT1, it held on
2026-08-02 but not on 2026-08-06 (rs12740374 ranks 2nd, ~15% below the top
variant, rs464218). This is expected to be fragile rather than a pipeline
defect: FOLD_0 excludes the SORT1 locus from training, and SORT1-specific
FOLD_0 predictions are independently documented in the manuscript (Fig. S5)
as sitting near the noise floor for this locus (RNA correlation "fell from
0.79 to -0.16"), so among 111 mostly-small, closely-spaced FOLD_0 scores,
rs12740374's lead over the runner-up was already narrow. `reproduction/
report.py::compare_figs2c` reports the exact per-gene rank for both the
archive and the fresh run rather than silently forcing a pass; it gates
`pass` on the more robust claim (top-3, not exact rank-1) plus a loose
correlation tolerance.

**Figure 1C itself** (main text, `ALL_FOLDS`) was not re-verified against
its original tight comparator tolerance (`reproduction/report.py::
compare_fig1c`, `rtol=1e-4, atol=1e-5`) under the 2026-08-06 backend state;
given the drift above, a fresh run today would likely not pass that
tolerance either, though the qualitative result (rs12740374 as the
dominant-effect variant) is unaffected at this magnitude of drift. This repo
does not currently loosen that comparator -- it is a main-text panel and
changing its pass/fail bar is a manuscript-level decision, not a
pipeline-maintenance one. A dedicated, deliberate full-repo freshness pass
(all figures, not just this one) closer to submission would give a
complete current-day picture.
