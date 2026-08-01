# Release checklist

## Panel identity

- [x] Freeze the author-approved manuscript text, supplementary tables and
      assembled SVGs with SHA-256 checksums (2026-08-01).
- [ ] Export the final manuscript PDF and supplementary PDF.
- [x] Trace every `provisional`/misattributed row's real generating script
      against the working archive (2026-08-01 audit; see MANIFEST_NOTES.md).
- [x] Resolve genuine content ambiguities (S1B, S7, S4B-D, S5A-D, S10, S11)
      with the author, rather than guessing (2026-08-01).
- [x] Reconcile figure numbering and composite assembly against the final
      manuscript package: four main figures and Figure S1–S10 (2026-08-01).
- [x] Preserve editable author-layout sources in the checksummed final
      composite SVGs under `figures/assembled/` and map each panel in the
      manifest (2026-08-01).
- [ ] Confirm that Figure 4 layout-only changes do not alter source tables or
      frozen statistics.

## Code and source data

- [ ] Copy/refactor only the analysis code used by selected panels
      (Figure 1C–F, Figure 2B/C/E/F, Figure 3E–G and S5B–C renderers are now
      ported; Figure 1B's API-backed compact export is complete; see
      `MANIFEST_NOTES.md`).
- [x] Write one compact source table per panel into `outputs/source_data/`
      (2026-08-01; 69 files, SHA-256 recorded in
      `outputs/run_manifests/source_data_sha256.tsv`). Exceptions: Fig1B (no
      cached table exists anywhere, needs a live API call) and large (>15 MB)
      tables recorded in `outputs/run_manifests/zenodo_pending_large_outputs.tsv`
      instead of being committed. Figure 1B was subsequently regenerated and
      is now included as a traced release-generated source directory.
- [ ] Add one deterministic rendering entry point per figure
      (`figures/fig1.py`, `fig2.py`, `fig3.py` and `figS5.py` now cover the
      source-backed panels described above; remaining panels retain exact
      legacy pointers).
- [x] Ensure figure rendering never requires an AlphaGenome API call
      (true for all ported panels; enforced by `figures/README.md` convention).
- [x] Record manuscript numbers in `outputs/manuscript_results.tsv`
      (257 source-verified rows generated from release tables on 2026-08-01).
- [ ] Run `python validate.py --require-release-assets` without errors
      (this strict command is the public-release gate; the default command is
      only the development-time structural check).

## Provenance

- [x] Resolve every source row in `data/SOURCES.tsv`, including exact study
      accessions/releases and honest non-redistribution records for bulk data.
- [x] Record URLs/accessions, release-audit access dates, retained-artifact
      checksums, genome builds, and terms.
- [x] Copy original frozen artifacts into `outputs/run_manifests/`
      (pre-existing from prior session).
- [x] Produce a SHA-256 manifest for compact AlphaGenome outputs
      (`outputs/run_manifests/source_data_sha256.tsv` +
      `outputs/run_manifests/zenodo_pending_large_outputs.tsv`, 2026-08-01).
- [x] Map all 69 compact release files to byte-identical source files in the
      frozen working archive (`source_data_provenance.tsv`, 2026-08-01).
- [x] Record coordinate conversions, allele orientation, score axes and
      length-preserving replacement conventions (`METHOD_CONVENTIONS.md`).

## Publication metadata

- [ ] Add the final manuscript title, authors, ORCIDs, and repository URL to
      `CITATION.cff`.
- [x] Confirm the code license and add `LICENSE` (MIT).
- [ ] Deposit large derived outputs on Zenodo and add the DOI to `README.md`.
- [ ] Replace the manuscript data-availability statement with the final GitHub
      and Zenodo locations.
- [ ] Confirm that redistribution of every included third-party table is
      permitted; otherwise include a fetch recipe rather than the file.

## Final validation

- [ ] Create the environment from `environment.yml` in a clean directory.
- [ ] Regenerate every figure from compact source tables.
- [ ] Compare regenerated SVG/PDF checksums or image renders with submission
      panels and document any harmless metadata-only differences.
- [ ] Run all numerical consistency checks and scan the manuscript for stale
      panel numbers, model names, tracks, and point estimates.
