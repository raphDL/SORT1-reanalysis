# SORT1 clean-room reproducibility report

- Run status: **complete**
- Reference comparison: **PASS**
- Started (UTC): `2026-08-02T19:56:48.916296+00:00`
- Repository commit: `f0423a73ae04fcd50b1cf915ab0d8fa5791e709f`
- Repository dirty during run: `True`
- Run directory: `/Volumes/T7/alphaGenome/repro_crash_test/runs/figure2_public_inputs_20260802T192638Z`
- AlphaGenome credential available during run: `True` (value never recorded)
- Timed step total: `7.2 s`
- Timed step total across retained attempts: `1256.9 s`
- Peak RSS reported by process: `412696576 bytes`
- Run-directory size: `4029987423 bytes`

## Scope

The analysis stage is isolated from `outputs/source_data/`. Published source tables are opened only after a successful analysis, by the comparison command.

## Post-run comparison code

- Path: `reproduction/report.py`
- Bytes: `18994`
- SHA-256: `b30af8dd717142f82f1762eb42598a122d3795c0fef466da0d59b10744adc56d`

## Executed reproduction code

| Path | Bytes | SHA-256 |
|---|---:|---|
| `reproduce.py` | 10565 | `39b6801ac4e2aa6cca892d3c5992abd7776146846eec5cc6507706cb2ff524f8` |
| `reproduction/__init__.py` | 72 | `ed9010e9708a159da43e98f7993611d2bf5ef32161c2f3c5381f6ee294ddd199` |
| `reproduction/common.py` | 13098 | `170f6eb17afd40a31538064de15d734c14aa44065978a5f374b5c34811f3444b` |
| `reproduction/figure1.py` | 24264 | `33bbcd2273448bbf62829ddb20c5e8addd0552795f4a22ff27b4d1cbb2d9117a` |
| `reproduction/figure2.py` | 43967 | `9b7bfd442e0bbd6d8cd0ff74e76f49704033281abb22db4d4a4197bf2f8eccde` |
| `reproduction/report.py` | 18994 | `b30af8dd717142f82f1762eb42598a122d3795c0fef466da0d59b10744adc56d` |

## Panel results

| Panel | Inputs | AlphaGenome regime | Scored units | API requests | Comparison |
|---|---|---|---:|---:|---|
| 2B | Wang 2018 spreadsheet + hg38 + AlphaGenome API | ALL_FOLDS | 52 | 7 | PASS |
| 2C | hg38 + AlphaGenome API | ALL_FOLDS | 98 | 13 | PASS |
| 2E | Kircher 2019 MPRA + AlphaGenome API | ALL_FOLDS | 5400 | 180 | PASS |
| 2F | Kircher 2019 MPRA + AlphaGenome API | ALL_FOLDS | 13335 | 447 | PASS |

## Timings

| Step | Status | Seconds |
|---|---|---:|
| Figure 2: download and checksum original public inputs | complete | 4.378 |
| 2B: reconstruct Wang repair products and score liver RNA | complete | 0.971 |
| 2C: construct and score the 11-by-11 deletion grid | complete | 0.334 |
| 2E: score the 1,798-substitution SORT1 Kircher benchmark | complete | 0.711 |
| 2F: score the six-element Kircher accessibility benchmark | complete | 0.829 |

## Prior attempts retained during resume

| Started UTC | Status | Timed seconds | Archive |
|---|---|---:|---|
| 2026-08-02T19:26:53.627234+00:00 | complete | 153.435 | `audit/attempts/2026-08-02T192653.627234_0000/` |
| 2026-08-02T19:29:47.582609+00:00 | failed | 90.632 | `audit/attempts/2026-08-02T192947.582609_0000/` |
| 2026-08-02T19:32:35.164301+00:00 | complete | 1005.607 | `audit/attempts/2026-08-02T193235.164301_0000/` |

## Downloads

| URL | Bytes | SHA-256 | Reused |
|---|---:|---|---|
| https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | 872949833 | `fb4243ebb014caf27111f24dd62b7ce42160f28581da6f8fcd6cba5977778d02` | True |
| manual_original_file:atv310103_ds.xls | 126976 | `99a2d4218e3e3a9afdd04ecded3a38878c001272f42eaadc99b8c3caeeefa049` | True |
| https://raw.githubusercontent.com/kircherlab/MPRA_SaturationMutagenesis/05d2ffb965090d3f5dd27dfb038cec493a15ab35/data/elements.tsv.gz | 2677485 | `fec2eed91fe27af3aae07ebce2eca65e9bad4bb6abba5d8c27f478887dd7b134` | True |

## Generated-file manifest

(Omitted here: 799 lines listing every generated file with its SHA-256 -- see run.json for the full audit trail. Available in the actual run directory this report was generated from.)

## Interpretation

A PASS means identities and counts were exact and numerical outputs met the explicit panel-specific equivalence thresholds recorded in `audit/comparison.json`. Figure 2F and non-live comparisons retain `rtol=1e-5, atol=1e-6`; live sequence/ISM panels additionally allow bounded sub-panel-unit API drift while requiring very high reference correlation. The observed maximum differences and plotted-correlation changes remain fully reported. The frozen tables were not used to generate any result.
