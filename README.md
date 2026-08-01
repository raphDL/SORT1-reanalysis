# SORT1 1p13.3 sequence-to-function reanalysis

This repository contains the analysis code, compact source data, run manifests,
and figure-generation code accompanying a computational reanalysis of the
`SORT1` 1p13.3 cholesterol locus with AlphaGenome. It reproduces the
sequence-to-function evidence chain — eQTL/tagging prioritization, observed
vs. predicted 3D contact, in-silico saturation mutagenesis, and cross-context
module-transfer/portability tests — reported in the accompanying manuscript.

## Repository status

The submission snapshot is frozen under `manuscript/` and
`figures/assembled/`, with four main figures and ten supplementary figures.
Supplementary numbering is final as Figure S1–S10 and is recorded in
`MANIFEST.tsv`. `manuscript/SUBMISSION_SNAPSHOT_SHA256.tsv` provides the
byte-level checksum manifest for the frozen text, tables and assembled SVGs.

Compact source data (with SHA-256 checksums) has been copied into
`outputs/source_data/` for every selected panel with a cached compact output;
the only remaining exceptions are four Zenodo-pending large tables. Figure 1B
was regenerated once from the authorized API and now has frozen numerical
arrays plus track-selection metadata.
Deterministic rendering scripts (`figures/fig*.py`) are ported for Figure 1C–F,
Figure 2B/C/E/F, Figure 3E–G, and the split Figure S5 panels B–C. These scripts
were test-rendered from the frozen compact tables without network or API
access. Figure 1B was recovered through a one-time authorized numerical export
because its legacy script did not persist the track arrays. See
`MANIFEST_NOTES.md` for the exact legacy generator for every panel. The working
archive is intentionally not copied wholesale—only files needed to reproduce
a reported result or render a published panel enter this repository.

## Layout

- `data/` — provenance records (`SOURCES.tsv`), download instructions, and
  small curated input tables that cannot be fetched automatically.
- `analysis/` — scripts that generate compact panel-level source tables from
  public inputs or AlphaGenome predictions (not all panels have one yet).
- `figures/` — deterministic rendering entry points as they are ported,
  author-approved editable composites in `figures/assembled/`, and
  checksummed canonical panel exports in `figures/rendered/`. The intended
  final interface is one script per manuscript figure; the current porting
  status is listed in `figures/README.md`.
- `outputs/source_data/` — one compact table (or small directory of tables)
  per panel, matching `MANIFEST.tsv`'s `release_source_table` column.
- `outputs/run_manifests/` — frozen decisions, model assignments, prespecified
  plans, and SHA-256 manifests, including
  `source_data_sha256.tsv` (checksums for everything in `outputs/source_data/`)
  and `zenodo_pending_large_outputs.tsv` (checksums for large derived tables
  that are documented but not committed to git; see below).
- `manuscript/` — frozen manuscript text, supplementary tables, snapshot
  metadata and SHA-256 checksums.
- `outputs/manuscript_results.tsv` — single source of truth for numbers
  quoted in the manuscript text (analysis ID + source-table pointer per row).
- `AUDIT_REPORT.md` — current release readiness, scientific corrections found
  during migration, and the remaining blockers to a public tag.

Root-level files provide the panel crosswalk (`MANIFEST.tsv` +
`MANIFEST_NOTES.md`), the frozen-analysis index (`PRESPECIFIED.md`), the
coordinate/allele/score conventions (`METHOD_CONVENTIONS.md`), the environment
spec, citation metadata, license, and the validation entry point.

## Installation

```
conda env create -f environment.yml
conda activate sort1-reanalysis
```

## Validating the release

```
python validate.py
```

Checks: required files present, manifest well-formed (no duplicate panels,
valid status values, valid AlphaGenome model-regime tokens), every file under
`outputs/source_data/` matches its recorded SHA-256, every Zenodo-pending
large output has a real (non-placeholder) SHA-256, no script under
`analysis/`/`figures/`/`src/` contains a personal absolute path, author-layout
(manual/adapted) panels carry composition notes, each compact table maps to a
byte-identical frozen-archive source, every AlphaGenome panel records its
model/input geometry, and the manuscript-results ledger points back to extant
release tables. Exits non-zero on structural or numerical-provenance defects.

Add `--workspace-root /path/to/alphaGenome` to additionally verify every
legacy-archive path still resolves (only meaningful if you have the full
38 GB working archive available; not needed to use or audit this release
repository on its own). Add `--require-release-assets` to additionally fail
on any selected panel still missing its rendered asset or source table. This
strict command, not the default structural check, is the public-release gate.

## Regenerating a figure

```
python figures/fig1.py            # renders Fig. 1C–F
python figures/fig2.py            # renders Fig. 2B/C/E/F
python figures/fig3.py            # renders the currently-ported Fig3 panels (E, F, G)
python figures/fig3.py --panel G  # renders a single panel
```

Figure scripts are pure, deterministic plotting from the compact tables in
`outputs/source_data/` — no AlphaGenome API calls, no public-data downloads,
no network access. Panels not yet ported (see above) must currently be
rendered from their legacy script in the working archive.

## Reproducibility model

The release separates two stages:

1. **Analysis** — potentially expensive AlphaGenome/API or public-data
   processing that writes a compact, versioned source table into
   `outputs/source_data/`.
2. **Rendering** — deterministic plotting from that compact source table,
   requiring no API access (`figures/fig*.py`).

Each published panel is mapped to its analysis script, rendering script,
source table, model regime, external datasets, and working-archive source in
`MANIFEST.tsv`.

## AlphaGenome

AlphaGenome predictions require authorized API access and are not automatically
executed by the figure-rendering commands. Model regime (`ALL_FOLDS`, `FOLD_0`,
or matched held-out fold) is recorded per panel in `MANIFEST.tsv`, and the
locus/model configuration is pinned in `config.yaml` and, because context
length differs by analysis, in
`outputs/run_manifests/panel_analysis_parameters.tsv`. The development
checkout used for the current analysis was:

- `alphagenome`: `327ea82371197812b337aed8e9e75df63ce5e429`
- `alphagenome_research`: `bb6a5276a51199d9589fc8929be6a9d946b7250b`

Figure 1B locus RNA-seq/CHIP-TF/DNase/ATAC arrays were regenerated once using
the ALL_FOLDS model and are now cached with their selected-track metadata and
hashes. The refactored exporter reads credentials only from the environment;
the hardcoded fallback discovered in the legacy working script was not copied.

## Data and code availability

Public-data provenance (accessions, URLs, genome builds, license/redistribution
terms and retained-artifact checksums) is tracked in `data/SOURCES.tsv`.
Large derived AlphaGenome-output tables that are
too large to commit (>15 MB; see `outputs/run_manifests/zenodo_pending_large_outputs.tsv`
for the current list, sizes, SHA-256 digests, and regenerating scripts) will
be deposited on Zenodo and linked by DOI before submission — **DOI: pending**.

## Citation and license

Author citation metadata is recorded in `CITATION.cff`; the final data DOI
will be added after the Zenodo deposit. Code is released under the license in
`LICENSE`; third-party datasets remain subject to their original terms as
recorded in `data/SOURCES.tsv`.

## A note on manually assembled panels

Panels marked `author-layout` in `MANIFEST.tsv` (schematics and adapted
graphical abstracts, e.g. Fig1A, Fig2A/D, Fig3D, Fig4A/D/I) were composed by
hand in Inkscape and are not reproducible by a Python plotting command. Their
editable author sources are preserved within the checksummed final composite
SVGs under `figures/assembled/`; the manifest identifies the relevant
composite for each such panel and records the source citation where needed.
