# SORT1 clean-room reproducibility report

- Run status: **complete**
- Reference comparison: **PASS**
- Started (UTC): `2026-08-05T12:40:48.615061+00:00`
- Repository commit: `81bb2bff81badee76053740e0673bd35b86a47ac`
- Repository dirty during run: `True`
- Run directory: `/Volumes/T7/alphaGenome/repro_crash_test/runs/figure4_bc_20260805`
- AlphaGenome credential available during run: `True` (value never recorded)
- Timed step total: `3760.5 s`
- Timed step total across retained attempts: `4030.3 s`
- Peak RSS reported by process: `5743624192 bytes`
- Run-directory size: `4636106304 bytes`

## Scope

The analysis stage is isolated from `outputs/source_data/`. Published source tables are opened only after a successful analysis, by the comparison command.

## Post-run comparison code

- Path: `reproduction/report.py`
- Bytes: `34166`
- SHA-256: `f7fb149a470fc43cd8b4d6b4f63d8c32ba92e9524c872543ee573c51ece6bc20`

## Executed reproduction code

| Path | Bytes | SHA-256 |
|---|---:|---|
| `reproduce.py` | 16005 | `089064815fc4ddd4d8a2c403c89e3c62fcae6b5ff251c4587f9636be99a3fbc1` |
| `reproduction/__init__.py` | 72 | `ed9010e9708a159da43e98f7993611d2bf5ef32161c2f3c5381f6ee294ddd199` |
| `reproduction/common.py` | 13098 | `170f6eb17afd40a31538064de15d734c14aa44065978a5f374b5c34811f3444b` |
| `reproduction/figure1.py` | 24311 | `08b93d844eaddd5524ae9be49862f7fff4b668839a8a9dde64bf59e8e2119220` |
| `reproduction/figure1_public.py` | 12261 | `30f23301e2d63630a888eca59847768b475c8ca004362e3a2a85afe25303d688` |
| `reproduction/figure2.py` | 45544 | `fb8e518cfcd4de6ecef4b6965f0a9ef772efb0b7ef4e6593a35df74cb7f51d8f` |
| `reproduction/figure3.py` | 72159 | `449c7ccc6ede6cd6c84b1cb1c42ff149786b7678a964e442f626920d53064c6c` |
| `reproduction/figure4.py` | 31939 | `0dddfa0dbdf6b0a94a59197b1848fc10b54b76241377e6631685fadc1db24e22` |
| `reproduction/report.py` | 33953 | `a300887beefe436246881dd93b8d6141f3f8e33f908b960157f729af113d7332` |

## Panel results

| Panel | Inputs | AlphaGenome regime | Scored units | API requests | Comparison |
|---|---|---|---:|---:|---|
| 4B | GRCh38 + frozen bottom100 recipient design + AlphaGenome API | ALL_FOLDS | 2500 | 20 | PASS |
| 4C | GRCh38 + HPA v24.1 + GENCODE v46 + AlphaGenome API | ALL_FOLDS | 5884 | 46 | PASS |

## Timings

| Step | Status | Seconds |
|---|---|---:|
| 4B: bottom100 315bp donor eight-distance sweep | complete | 591.591 |
| 4B/4C: download UCSC-sourced hg38 (the build the working archive actually used) | complete | 4.381 |
| 4B: build 315bp donor + bottom100 distance-sweep design | complete | 2.036 |
| 4B: score states 1-128 of 2500 | complete | 37.304 |
| 4B: score states 129-256 of 2500 | complete | 25.235 |
| 4B: score states 257-384 of 2500 | complete | 26.772 |
| 4B: score states 385-512 of 2500 | complete | 28.851 |
| 4B: score states 513-640 of 2500 | complete | 24.919 |
| 4B: score states 641-768 of 2500 | complete | 30.001 |
| 4B: score states 769-896 of 2500 | complete | 31.006 |
| 4B: score states 897-1024 of 2500 | complete | 28.794 |
| 4B: score states 1025-1152 of 2500 | complete | 29.019 |
| 4B: score states 1153-1280 of 2500 | complete | 31.501 |
| 4B: score states 1281-1408 of 2500 | complete | 25.957 |
| 4B: score states 1409-1536 of 2500 | complete | 23.449 |
| 4B: score states 1537-1664 of 2500 | complete | 26.455 |
| 4B: score states 1665-1792 of 2500 | complete | 26.700 |
| 4B: score states 1793-1920 of 2500 | complete | 24.127 |
| 4B: score states 1921-2048 of 2500 | complete | 24.051 |
| 4B: score states 2049-2176 of 2500 | complete | 24.162 |
| 4B: score states 2177-2304 of 2500 | complete | 63.874 |
| 4B: score states 2305-2432 of 2500 | complete | 28.166 |
| 4B: score states 2433-2500 of 2500 | complete | 14.666 |
| 4B: compute deltas and summarize by distance | complete | 0.211 |
| 4C: HPA bottom/middle/top-500 cohort fold-change at 30bp | complete | 1310.803 |
| 4B/4C: download UCSC-sourced hg38 (the build the working archive actually used) | complete | 4.330 |
| 4C: stage manually supplied HPA v24.1 consensus tissue RNA expression | complete | 0.216 |
| 4C: download GENCODE v46 annotation feather | complete | 16.192 |
| 4C: select HPA bottom/middle/top-500 cohorts and build 30bp donor-transfer design | complete | 8.448 |
| 4C: score states 1-128 of 5884 | complete | 30.939 |
| 4C: score states 129-256 of 5884 | complete | 24.986 |
| 4C: score states 257-384 of 5884 | complete | 25.144 |
| 4C: score states 385-512 of 5884 | complete | 25.646 |
| 4C: score states 513-640 of 5884 | complete | 25.146 |
| 4C: score states 641-768 of 5884 | complete | 24.973 |
| 4C: score states 769-896 of 5884 | complete | 25.147 |
| 4C: score states 897-1024 of 5884 | complete | 25.346 |
| 4C: score states 1025-1152 of 5884 | complete | 26.405 |
| 4C: score states 1153-1280 of 5884 | complete | 24.473 |
| 4C: score states 1281-1408 of 5884 | complete | 24.702 |
| 4C: score states 1409-1536 of 5884 | complete | 26.375 |
| 4C: score states 1537-1664 of 5884 | complete | 25.731 |
| 4C: score states 1665-1792 of 5884 | complete | 24.926 |
| 4C: score states 1793-1920 of 5884 | complete | 25.984 |
| 4C: score states 1921-2048 of 5884 | complete | 24.698 |
| 4C: score states 2049-2176 of 5884 | complete | 26.273 |
| 4C: score states 2177-2304 of 5884 | complete | 25.935 |
| 4C: score states 2305-2432 of 5884 | complete | 25.280 |
| 4C: score states 2433-2560 of 5884 | complete | 26.752 |
| 4C: score states 2561-2688 of 5884 | complete | 32.356 |
| 4C: score states 2689-2816 of 5884 | complete | 26.454 |
| 4C: score states 2817-2944 of 5884 | complete | 25.861 |
| 4C: score states 2945-3072 of 5884 | complete | 41.824 |
| 4C: score states 3073-3200 of 5884 | complete | 26.797 |
| 4C: score states 3201-3328 of 5884 | complete | 25.496 |
| 4C: score states 3329-3456 of 5884 | complete | 25.742 |
| 4C: score states 3457-3584 of 5884 | complete | 26.110 |
| 4C: score states 3585-3712 of 5884 | complete | 26.232 |
| 4C: score states 3713-3840 of 5884 | complete | 25.719 |
| 4C: score states 3841-3968 of 5884 | complete | 28.574 |
| 4C: score states 3969-4096 of 5884 | complete | 25.978 |
| 4C: score states 4097-4224 of 5884 | complete | 27.239 |
| 4C: score states 4225-4352 of 5884 | complete | 26.397 |
| 4C: score states 4353-4480 of 5884 | complete | 25.659 |
| 4C: score states 4481-4608 of 5884 | complete | 26.743 |
| 4C: score states 4609-4736 of 5884 | complete | 29.891 |
| 4C: score states 4737-4864 of 5884 | complete | 29.916 |
| 4C: score states 4865-4992 of 5884 | complete | 34.925 |
| 4C: score states 4993-5120 of 5884 | complete | 31.163 |
| 4C: score states 5121-5248 of 5884 | complete | 26.731 |
| 4C: score states 5249-5376 of 5884 | complete | 27.697 |
| 4C: score states 5377-5504 of 5884 | complete | 27.354 |
| 4C: score states 5505-5632 of 5884 | complete | 25.954 |
| 4C: score states 5633-5760 of 5884 | complete | 26.603 |
| 4C: score states 5761-5884 of 5884 | complete | 28.659 |
| 4C: compute fold changes by cohort | complete | 0.325 |

## Prior attempts retained during resume

| Started UTC | Status | Timed seconds | Archive |
|---|---|---:|---|
| 2026-08-05T12:35:46.871591+00:00 | running | 269.812 | `audit/attempts/2026-08-05T123546.871591_0000/` |

## Downloads

| URL | Bytes | SHA-256 | Reused |
|---|---:|---|---|
| https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz | 983659424 | `c1dd87068c254eb53d944f71e51d1311964fce8de24d6fc0effc9c61c01527d4` | True |
| supplied:<local-machine>/alphaGenome/investigation/SORT1_cholesterol_musunuru_2010/report/panel_scramble_no_expression/HPA_rna_tissue_consensus.tsv | 39004810 | `cdedaeaf3cdfc89e22b3891ea24ae2afabc0afd26d8883076121a363608450b6` | False |
| https://storage.googleapis.com/alphagenome/reference/gencode/hg38/gencode.v46.annotation.gtf.gz.feather | 333040258 | `7b10f643d96e1142ef058d9c08487f4360552cea29c8339f195a5a96489dbb4c` | False |

## Generated-file manifest

(Omitted here: 8,402 lines listing every generated file with its SHA-256 -- see run.json for the full audit trail. Available in the actual run directory this report was generated from.)

## Interpretation

A PASS means identities and counts were exact and numerical outputs met the explicit panel-specific equivalence thresholds recorded in `audit/comparison.json`. Figure 2F and non-live comparisons retain `rtol=1e-5, atol=1e-6`; live sequence/ISM panels additionally allow bounded sub-panel-unit API drift while requiring very high reference correlation. The observed maximum differences and plotted-correlation changes remain fully reported. The frozen tables were not used to generate any result.
