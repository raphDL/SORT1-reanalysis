# SORT1 clean-room reproducibility report

> **Note**: this run was intentionally stopped by user decision after 256 of ~900,000 planned predictions (see REPRODUCIBILITY_NEXT_STEPS.md, "Figure 4H" update) -- not a failure or crash. Included here for transparency about exactly how far the run got before stopping.


- Run status: **running**
- Reference comparison: **NOT YET PASSED**
- Started (UTC): `2026-08-05T21:16:44.888286+00:00`
- Repository commit: `99e99318f06bb4a56f67168589ac488ca83d04bf`
- Repository dirty during run: `True`
- Run directory: `/Volumes/T7/alphaGenome/repro_crash_test/runs/figure4h_20260805`
- AlphaGenome credential available during run: `True` (value never recorded)
- Timed step total: `203.9 s`
- Timed step total across retained attempts: `203.9 s`
- Peak RSS reported by process: `1612005376 bytes`
- Run-directory size: `4257244484 bytes`

## Scope

The analysis stage is isolated from `outputs/source_data/`. Published source tables are opened only after a successful analysis, by the comparison command.

## Executed reproduction code

| Path | Bytes | SHA-256 |
|---|---:|---|
| `reproduce.py` | 17231 | `dbfa4dba6b8e50115c6a523a84044c25c7d9de8557d69c6b8f7513ba91251158` |
| `reproduction/__init__.py` | 72 | `ed9010e9708a159da43e98f7993611d2bf5ef32161c2f3c5381f6ee294ddd199` |
| `reproduction/common.py` | 13098 | `170f6eb17afd40a31538064de15d734c14aa44065978a5f374b5c34811f3444b` |
| `reproduction/figure1.py` | 24311 | `08b93d844eaddd5524ae9be49862f7fff4b668839a8a9dde64bf59e8e2119220` |
| `reproduction/figure1_public.py` | 12261 | `30f23301e2d63630a888eca59847768b475c8ca004362e3a2a85afe25303d688` |
| `reproduction/figure2.py` | 45544 | `fb8e518cfcd4de6ecef4b6965f0a9ef772efb0b7ef4e6593a35df74cb7f51d8f` |
| `reproduction/figure3.py` | 72159 | `449c7ccc6ede6cd6c84b1cb1c42ff149786b7678a964e442f626920d53064c6c` |
| `reproduction/figure4.py` | 31939 | `0dddfa0dbdf6b0a94a59197b1848fc10b54b76241377e6631685fadc1db24e22` |
| `reproduction/figure4ef.py` | 48018 | `ec152c36518d284b143e9dd987054b4e7a0b0c929c806299f4c822df61cbfe27` |
| `reproduction/figure4g.py` | 4986 | `a09beffeb4fb4fe64537c5a3149a761c04eb22948f568a61f15d0f36452e3d62` |
| `reproduction/figure4h.py` | 12187 | `672fd635827a408600406f9f8bd9d0e8a2d762826d241840cf1e230ac0397397` |
| `reproduction/report.py` | 42724 | `8a7ad7487879cc376e233c1cd15e0bfabdaf85ffd714b0e27ef4b8716ffed06c` |

## Panel results

| Panel | Inputs | AlphaGenome regime | Scored units | API requests | Comparison |
|---|---|---|---:|---:|---|
| 4H | GRCh38 + AlphaGenome API (exhaustive +/-50kb x 3-alt x 3-tissue ISM) | ALL_FOLDS | 256 | 8 | NOT COMPARED |

## Timings

| Step | Status | Seconds |
|---|---|---:|
| 4H: exhaustive +/-50kb x 3-alt x 3-tissue regional ISM synergy scan | running | 0.000 |
| 4B/4C: download UCSC-sourced hg38 (the build the working archive actually used) | complete | 124.432 |
| 4H: build the +/-50kb x 3-alt candidate-variant universe | complete | 0.029 |
| 4H liver: score variants 1-32 of 300003 (of 300003 total) | complete | 9.991 |
| 4H liver: score variants 33-64 of 300003 (of 300003 total) | complete | 9.508 |
| 4H liver: score variants 65-96 of 300003 (of 300003 total) | complete | 9.365 |
| 4H liver: score variants 97-128 of 300003 (of 300003 total) | complete | 9.757 |
| 4H liver: score variants 129-160 of 300003 (of 300003 total) | complete | 9.735 |
| 4H liver: score variants 161-192 of 300003 (of 300003 total) | complete | 12.308 |
| 4H liver: score variants 193-224 of 300003 (of 300003 total) | complete | 9.228 |
| 4H liver: score variants 225-256 of 300003 (of 300003 total) | complete | 9.521 |
| 4H liver: score variants 257-288 of 300003 (of 300003 total) | running | 0.000 |

## Downloads

| URL | Bytes | SHA-256 | Reused |
|---|---:|---|---|
| https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz | 983659424 | `c1dd87068c254eb53d944f71e51d1311964fce8de24d6fc0effc9c61c01527d4` | False |

## Generated-file manifest

(Omitted here: 9 lines listing every generated file with its SHA-256.)

## Interpretation

A PASS means identities and counts were exact and numerical outputs met the explicit panel-specific equivalence thresholds recorded in `audit/comparison.json`. Figure 2F and non-live comparisons retain `rtol=1e-5, atol=1e-6`; live sequence/ISM panels additionally allow bounded sub-panel-unit API drift while requiring very high reference correlation. The observed maximum differences and plotted-correlation changes remain fully reported. The frozen tables were not used to generate any result.
