# SORT1 clean-room reproducibility report

- Run status: **complete**
- Reference comparison: **PASS**
- Started (UTC): `2026-08-05T11:26:31.664464+00:00`
- Repository commit: `23171f0e1f10e92f740c7a15dbb635e4876a7fdb`
- Repository dirty during run: `True`
- Run directory: `/Volumes/T7/alphaGenome/repro_crash_test/runs/figure3_all_panels_20260805`
- AlphaGenome credential available during run: `True` (value never recorded)
- Timed step total: `434.4 s`
- Timed step total across retained attempts: `434.4 s`
- Peak RSS reported by process: `3398008832 bytes`
- Run-directory size: `4091507392 bytes`

## Scope

The analysis stage is isolated from `outputs/source_data/`. Published source tables are opened only after a successful analysis, by the comparison command.

## Post-run comparison code

- Path: `reproduction/report.py`
- Bytes: `31876`
- SHA-256: `1aac1febeb9fd558d1323a52e955b8999fabb195c461556948d38fa2f34ff9e3`

## Executed reproduction code

| Path | Bytes | SHA-256 |
|---|---:|---|
| `reproduce.py` | 15085 | `bf54bc4f48568b9eac13c6ff385792f40bc27c9fce5ae60336bc3d18b5766b18` |
| `reproduction/__init__.py` | 72 | `ed9010e9708a159da43e98f7993611d2bf5ef32161c2f3c5381f6ee294ddd199` |
| `reproduction/common.py` | 13098 | `170f6eb17afd40a31538064de15d734c14aa44065978a5f374b5c34811f3444b` |
| `reproduction/figure1.py` | 24311 | `08b93d844eaddd5524ae9be49862f7fff4b668839a8a9dde64bf59e8e2119220` |
| `reproduction/figure1_public.py` | 12261 | `30f23301e2d63630a888eca59847768b475c8ca004362e3a2a85afe25303d688` |
| `reproduction/figure2.py` | 45544 | `fb8e518cfcd4de6ecef4b6965f0a9ef772efb0b7ef4e6593a35df74cb7f51d8f` |
| `reproduction/figure3.py` | 72159 | `449c7ccc6ede6cd6c84b1cb1c42ff149786b7678a964e442f626920d53064c6c` |
| `reproduction/report.py` | 31876 | `1aac1febeb9fd558d1323a52e955b8999fabb195c461556948d38fa2f34ff9e3` |

## Panel results

| Panel | Inputs | AlphaGenome regime | Scored units | API requests | Comparison |
|---|---|---|---:|---:|---|
| 3A | GRCh38 + AlphaGenome API | ALL_FOLDS | 0 | 0 | PASS |
| 3B | GRCh38 + AlphaGenome API | ALL_FOLDS | 0 | 0 | PASS |
| 3C | Figure 3B outputs + JASPAR 2024 CORE + GRCh38 | none (local PWM scan) | 0 | 0 | PASS |
| 3E | GRCh38 + AlphaGenome API | ALL_FOLDS | 0 | 0 | PASS |
| 3F | GRCh38 + AlphaGenome API | ALL_FOLDS | 0 | 0 | PASS |
| 3G | GRCh38 + AlphaGenome API | ALL_FOLDS | 0 | 0 | PASS |

## Timings

| Step | Status | Seconds |
|---|---|---:|
| Figure 3: download and checksum GRCh38 | complete | 4.572 |
| 3B: native-locus 501-bp three-gene ISM | complete | 1.479 |
| 3C: scan JASPAR motifs and score PWM compatibility vs Figure 3B ISM | complete | 45.311 |
| 3C: download JASPAR 2024 CORE PFM matrices | complete | 1.240 |
| 3C: scan JASPAR motifs and score PWM compatibility vs Figure 3B ISM | complete | 42.395 |
| 3G: expanded component-necessity audit | complete | 4.317 |
| 3G: build component-necessity scramble design | complete | 2.677 |
| 3G: compute retention and render | complete | 0.058 |
| 3E: directional single-arm motif-protected recovery | complete | 8.115 |
| 3E: build directional single-arm recovery design | complete | 6.260 |
| 3E: compute retention and render | complete | 0.065 |
| 3A: 100kb regional two-stage RNA(TSS) ISM scan | complete | 1.862 |
| 3A: write derived tables and render | complete | 0.437 |
| 3F: wide-main-panel 1bp boundary grid | complete | 170.886 |
| 3F: build wide-main-panel 1bp boundary grid design | complete | 143.905 |
| 3F: compute retention surfaces and render | complete | 0.862 |

## Downloads

| URL | Bytes | SHA-256 | Reused |
|---|---:|---|---|
| https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | 872949833 | `fb4243ebb014caf27111f24dd62b7ce42160f28581da6f8fcd6cba5977778d02` | True |
| https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt | 283814 | `2a6c7c24afb0614ed418c9f02b68845adebf38f80f6cebabbcd5de804eaacb59` | False |

## Generated-file manifest

(Omitted here: 17,460 lines listing every generated file with its SHA-256 -- see run.json for the full audit trail. Available in the actual run directory this report was generated from.)

## Interpretation

A PASS means identities and counts were exact and numerical outputs met the explicit panel-specific equivalence thresholds recorded in `audit/comparison.json`. Figure 2F and non-live comparisons retain `rtol=1e-5, atol=1e-6`; live sequence/ISM panels additionally allow bounded sub-panel-unit API drift while requiring very high reference correlation. The observed maximum differences and plotted-correlation changes remain fully reported. The frozen tables were not used to generate any result.
