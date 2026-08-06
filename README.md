# SORT1 1p13.3 sequence-to-function reanalysis

This repository contains the analysis code, compact source data, run manifests,
and figure-generation code accompanying a computational reanalysis of the
`SORT1` 1p13.3 cholesterol locus with AlphaGenome. It reproduces the
sequence-to-function evidence chain — eQTL/tagging prioritization, observed
vs. predicted 3D contact, in-silico saturation mutagenesis, and cross-context
module-transfer/portability tests — reported in the accompanying manuscript.

## Reproducibility capsule

**To reproduce the figures and data from scratch — no access to the
authors' working archive, no reused predictions — go to
[`reproduce.py`](reproduce.py) and [`reproduction/`](reproduction/).**
That is the entire clean-room capsule: a single CLI (`doctor` / `prepare` /
`run` / `compare`) that downloads every public input fresh, calls the
AlphaGenome API fresh for every prediction, and checksums the result
against this repository's own frozen publication tables in
`outputs/source_data/`. See "Reproducing Figure N analyses from scratch"
below for the exact commands per figure, and the "Reproducibility model"
section for what each of the four subcommands does and does not touch.

Evidence that this capsule has actually been run end-to-end — not just
written — is in [`reproduction_example_reruns/`](reproduction_example_reruns/):
audit trails (timestamps, exact commit, per-panel PASS/FAIL, timings, API
call counts) from real executions on a separate machine/volume from the
development checkout, each starting from an empty run directory so every
download and every prediction in them was freshly fetched or scored for
that run.

## Repository status

The submission snapshot is frozen under `manuscript/` and
`figures/assembled/`, with four main figures and ten supplementary figures.
Supplementary numbering is final as Figure S1–S10 and is recorded in
`MANIFEST.tsv`. `manuscript/SUBMISSION_SNAPSHOT_SHA256.tsv` provides the
byte-level checksum manifest for the frozen text, tables and assembled SVGs.

Frozen publication-reference data (with SHA-256 checksums) is stored under
`outputs/source_data/` for every selected panel with a cached compact output;
the only remaining exceptions are four Zenodo-pending large tables. These
tables support deterministic rendering and post-run comparison, but are not
inputs to the clean-room analysis runner.
Deterministic rendering scripts (`figures/fig*.py`) are ported for Figure 1C–F,
Figure 2B/C/E/F, Figure 3E–G, and the split Figure S5 panels B–C. These scripts
were test-rendered from the frozen compact tables without network or API
access. The initial Figure 1B reference was recovered through an authorized
numerical export because its legacy script did not persist the track arrays;
the clean-room runner can now regenerate it independently. See
`MANIFEST_NOTES.md` for the exact legacy generator for every panel. The working
archive is intentionally not copied wholesale—only files needed to reproduce
a reported result or render a published panel enter this repository.

## Layout

- `data/` — provenance records (`SOURCES.tsv`), download instructions, and
  small curated input tables that cannot be fetched automatically.
- `analysis/` — scripts that generate compact panel-level source tables from
  public inputs or AlphaGenome predictions (not all panels have one yet).
- `reproduction/` plus `reproduce.py` — isolated analysis, public-input
  download, timing, checksum, comparison, and audit-report machinery.
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

## Reproducing Figure 1 analyses from scratch

The clean-room entry point regenerates analytical outputs without reading the
committed tables under `outputs/source_data/`. For the first supported vertical
slice it reproduces Figure 1B, Figure 1C's middle column, and Figure 1E:

- **Figure 1B:** one AlphaGenome `ALL_FOLDS` locus-track prediction;
- **Figure 1C middle:** reconstruct the 111 variants from the full GLGC 2013
  download, the UCSC hg19-to-hg38 chain (from a checksum-identical public
  mirror), and the NCBI GRCh38 no-alt analysis-set FASTA, then run the
  recommended AlphaGenome RNA exon-mask scorer for every variant;
- **Figure 1E:** one held-out AlphaGenome `FOLD_0` contact-map prediction.

Create a credential file locally (it is ignored by git):

```bash
cp .env.example .env
# Edit .env and replace the placeholder with your authorized key.
python reproduce.py doctor
```

Alternatively, use the visible text-file template (useful in Finder):

```bash
cp ALPHAGENOME_API_KEY.example.txt ALPHAGENOME_API_KEY.txt
# Edit ALPHAGENOME_API_KEY.txt, keeping the ALPHAGENOME_API_KEY="..." format.
python reproduce.py doctor
```

`ALPHAGENOME_API_KEY.txt` is also ignored by Git and is discovered
automatically. Use either it or `.env`; neither credential value is printed or
written to an audit report.

Start each audit in a new, empty directory. A removable disk may be used for
downloads and run products while the Conda environment remains on an APFS disk:

```bash
# Optional: download and validate public inputs before an API key is available.
python reproduce.py prepare \
  --run-dir /Volumes/T7/alphaGenome/repro_crash_test/runs/figure1_clean

python reproduce.py run \
  --panels 1B,1C,1D,1E,1F \
  --run-dir /Volumes/T7/alphaGenome/repro_crash_test/runs/figure1_clean \
  --resume

python reproduce.py compare \
  --run-dir /Volumes/T7/alphaGenome/repro_crash_test/runs/figure1_clean
```

`run` refuses a non-empty directory unless `--resume` is supplied. Interrupted
Figure 1C batches are checkpointed. The run records input URLs and checksums,
package and Git versions, model regimes, request units, retries, timings, peak
RSS, disk use, and generated-file hashes in `audit/run.json` and
`audit/REPRODUCIBILITY_REPORT.md`. The API key itself is never recorded.

Only `compare` reads the frozen publication tables. It applies explicit numeric
tolerances and updates the report with a per-panel PASS/FAIL result. Supplying
`--max-variants` is available for development smoke tests, but such a partial
Figure 1C run cannot pass the publication comparison.

For a key file stored elsewhere, place the global option before the subcommand.
Both a raw key and an `ALPHAGENOME_API_KEY="..."` assignment are accepted:

```bash
python reproduce.py --api-key-file /secure/path/alphagenome.key doctor
```

This reconstructs Figure 1C-left from GTEx v7 liver cis-eQTLs, Figure 1C-right
from phased 1000 Genomes EUR haplotypes, and Figure 1D/1F from public
4DNFICSTCJQZ HepG2 Hi-C. The runner downloads GTEx (about 3.4 GB), uses
`bcftools` remote region access for the much larger chromosome-1 VCF, and uses
HTTP byte-range access for the 14.7-GB `.hic` file. Users who prefer explicit
downloads can supply the original files:

```bash
python reproduce.py run --panels 1B,1C,1D,1E,1F \
  --gtex-file /data/Liver.allpairs.txt.gz \
  --onekg-vcf /data/1kgp_chr1_sort1.vcf.gz \
  --onekg-panel /data/integrated_call_samples_v3.20130502.ALL.panel \
  --hic-file /data/4DNFICSTCJQZ.hic \
  --run-dir /data/figure1_clean
```

The supplied VCF may be the full phased chr1 release or a `bcftools view`
interval covering `chr1:109209432-109340504`; keep its `.tbi` alongside it.
Checksums, source URLs, byte counts, and whether remote range access was used
are recorded in the audit. Figure 1A remains an author-layout schematic.

## Reproducing Figure 2 analyses from scratch

Figure 2B, 2C, 2E, and 2F are supported by the same isolated runner. Figure
2A and 2D are literature schematics rather than computational panels.

## Reproducing Figure 3 analyses from scratch

Figure 3A, 3B, 3C, 3E, 3F, and 3G are all supported as checkpointed
clean-room AlphaGenome runs. Figure 3B constructs the intact rs12740374-T
sequence from downloaded GRCh38, predicts the intact baseline and every
three-way substitution across the 501-bp native locus, and derives the
three-gene positional loss profile and hotspots without reading
publication tables; 3C is a deterministic downstream JASPAR 2024 analysis
of the 3B substitutions (no additional API calls); 3A is a two-stage
regional (100kb) coordinated-RNA ISM scan; 3E/3F/3G are directional
single-arm motif-protected recovery, a wide boundary grid, and an expanded
component-necessity audit, respectively. Figure 3D is an author-layout
schematic.

```bash
python reproduce.py --api-key-file /secure/path/alphagenome.key run \
  --panels 3A,3B,3C,3E,3F,3G --batch-size 24 --max-workers 4 \
  --run-dir /data/figure3_clean
python reproduce.py compare --run-dir /data/figure3_clean
```

3C and 3A/3E/3F/3G each depend on 3B's native-locus ISM output; `run`
regenerates that prerequisite automatically if it isn't already present in
the run directory.

## Reproducing Figure 4 analyses from scratch

Figure 4B, 4C, 4E, 4F, and 4G are supported as checkpointed clean-room
AlphaGenome runs, all real and freshly scored:

- **4B:** a 315bp asymmetric donor (rs12740374 major/minor allele core)
  transferred to 100 low-expression liver recipient TSSs across eight
  upstream distances;
- **4C:** the same donor-transfer design at a fixed 30bp distance across
  HPA v24.1 liver-expression bottom/middle/top-500 gene cohorts (fully
  re-derived from the public HPA download at run time, not frozen);
- **4E/4F:** a chromosome-1, Hi-C-guided high/low distal-contact-site
  selection (from the public 4DN HepG2 Hi-C map, zero AlphaGenome cost)
  followed by real 315bp T/G transfer scoring at each selected site
  (13,517 predictions);
- **4G:** a single AlphaGenome variant-scoring call (rs12740374, RNA_SEQ)
  across tissue ontologies.

```bash
python reproduce.py --api-key-file /secure/path/alphagenome.key run \
  --panels 4B,4C,4E,4F,4G --batch-size 32 --max-workers 8 \
  --run-dir /data/figure4_clean
python reproduce.py compare --run-dir /data/figure4_clean
```

Figure 4H (an exhaustive +/-50kb x 3-allele x 3-tissue single-nucleotide
mutagenesis scan, ~900,000 predictions) is ported (`--panels 4H`) but its
full run has not been completed — see `REPRODUCIBILITY_NEXT_STEPS.md` for
why, including a real audited finding that individual-SNP effects at that
resolution sit at or below AlphaGenome's own measurement drift over time.
Its reference table, like Figure 4J's (a TF-motif-insertion discovery map,
~1.6M predicted variants, scoped but not yet ported), is Zenodo-pending
rather than committed to this repository; both accept a manually supplied,
checksum-verified copy (`compare --reference-4h-file ...`) in the
meantime. Figure 4A, 4D, and 4I are author-layout schematics.

- **2B:** reconstruct 50 unique repair products from the original Wang et al.
  spreadsheet, place them on the rs12740374-T background, and predict liver
  RNA for each complete edited sequence;
- **2C:** construct all 121 upstream-by-downstream deletion geometries from
  downloaded GRCh38 and predict each complete sequence;
- **2E:** download the Kircher et al. release and score 1,798 substitutions on
  the rs12740374-T construct background for HepG2 RNA, ATAC, and H3K27ac;
- **2F:** apply the prespecified barcode filter and score cell-line-matched
  DNase (six elements) and ATAC (four elements).

The Kircher file downloads automatically from an immutable upstream Git
commit. The Wang publisher currently places its spreadsheet download behind
bot protection. If automatic preparation reports that problem, download the
original `atv310103_ds.xls` from DOI `10.1161/ATVBAHA.117.310103`; the runner
accepts it only if its SHA-256 is
`99a2d4218e3e3a9afdd04ecded3a38878c001272f42eaadc99b8c3caeeefa049`.

```bash
python reproduce.py --api-key-file /secure/path/alphagenome.key prepare \
  --panels 2B,2C,2E,2F \
  --wang-xls /path/to/atv310103_ds.xls \
  --run-dir /Volumes/T7/alphaGenome/repro_crash_test/runs/figure2_clean

python reproduce.py --api-key-file /secure/path/alphagenome.key run \
  --panels 2B,2C,2E,2F \
  --wang-xls /path/to/atv310103_ds.xls \
  --run-dir /Volumes/T7/alphaGenome/repro_crash_test/runs/figure2_clean \
  --resume

python reproduce.py compare \
  --run-dir /Volumes/T7/alphaGenome/repro_crash_test/runs/figure2_clean
```

The audit distinguishes actual API requests from scored units (complete
sequences or variant/modality scores), checkpoints every batch, and reports
per-panel numerical PASS/FAIL against the frozen publication outputs only in
the final `compare` stage. The comparison requires exact sequence/variant
identities and row counts. It reports strict differences for every numerical
field and uses prespecified panel-unit equivalence bounds for live API outputs
(plus correlation requirements), because repeated remote inference can show
small floating-point drift. The exact thresholds and every observed maximum
difference are written to `audit/comparison.json`; they are never silently
rounded away.

## Rendering from frozen publication tables

```
python figures/fig1.py            # renders Fig. 1C–F
python figures/fig2.py            # renders Fig. 2B/C/E/F
python figures/fig3.py            # renders the currently-ported Fig3 panels (E, F, G)
python figures/fig3.py --panel G  # renders a single panel
```

These figure scripts are pure, deterministic plotting from the compact tables in
`outputs/source_data/` — no AlphaGenome API calls, no public-data downloads,
no network access. Panels not yet ported (see above) must currently be
rendered from their legacy script in the working archive.

## Reproducibility model

The release separates three stages:

1. **Clean-room analysis** — potentially expensive AlphaGenome/API or
   public-data processing that writes new predictions and derived tables into
   an isolated run directory. It never writes into or reads from the frozen
   publication-reference tables.
2. **Rendering** — deterministic plotting from that compact source table,
   requiring no API access.
3. **Comparison** — an explicitly separate post-run operation that may read
   `outputs/source_data/` to test the new numerical results against the frozen
   publication reference.

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

The frozen Figure 1B reference contains locus RNA-seq/CHIP-TF/DNase/ATAC arrays
and selected-track metadata. `reproduce.py` regenerates those arrays from the
API in a new run directory. Credentials are read from `.env`, the process
environment, or an explicitly supplied one-line key file; the hardcoded
fallback discovered in the legacy working script was not copied.

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
