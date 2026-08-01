# Figure rendering

Figure scripts render the published panels from compact files in
`outputs/source_data/`. They must not call AlphaGenome or download public data.

The intended public interface is one entry point per manuscript figure
(`fig1.py` through `fig4.py`, and `figS1.py` through `figS11.py`). Vector SVG is
the archival output; PDF and raster copies may be generated for submission or
visual comparison. Author-composed schematics are stored as editable SVGs and
are not regenerated from Python.

Currently ported deterministic entry points:

- `fig1.py`: Figure 1C, 1D, 1E and 1F;
- `fig2.py`: Figure 2B, 2C, 2E and 2F;
- `fig3.py`: Figure 3E, 3F and 3G;
- `figS5.py`: separately rendered ATAC (S5B) and H3K27ac (S5C) panels.

Figure 1B required a one-time authorized ALL_FOLDS regeneration because its
legacy script plotted numerical tracks directly from an in-memory AlphaGenome
response. The regenerated arrays and selected-track metadata are now retained
under `outputs/source_data/Figure1B_locus_tracks/`; the release contains no API
credential.

The Python entry points write SVG, PDF and PNG outputs to `figures/rendered/`
by default. Use `--output-dir` for a non-destructive comparison render. The
frozen submitted composites under `figures/assembled/` remain the authority
for author-positioned schematic panels and final multi-panel layout.

All other staged panel assets are byte-identical copies of the canonical
working-archive export and retain a provenance row in
`outputs/run_manifests/release_asset_provenance.tsv` until their plotting code
is ported.
