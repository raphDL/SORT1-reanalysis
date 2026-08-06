# SORT1 clean-room reproducibility report

- Run status: **complete**
- Reference comparison: **PASS**
- Started (UTC): `2026-08-05T18:24:34.874411+00:00`
- Repository commit: `328cd24b866deb2825fdbf8115e359f6f08cb9ce`
- Repository dirty during run: `True`
- Run directory: `/Volumes/T7/alphaGenome/repro_crash_test/runs/figure4g_20260805`
- AlphaGenome credential available during run: `True` (value never recorded)
- Timed step total: `1.5 s`
- Timed step total across retained attempts: `1.5 s`
- Peak RSS reported by process: `229474304 bytes`
- Run-directory size: `6880 bytes`

## Scope

The analysis stage is isolated from `outputs/source_data/`. Published source tables are opened only after a successful analysis, by the comparison command.

## Post-run comparison code

- Path: `reproduction/report.py`
- Bytes: `40063`
- SHA-256: `43159c9aac0afd3a70c46d5c438b4e6b4f22b952b34d2df03918ea534ce94c43`

## Executed reproduction code

| Path | Bytes | SHA-256 |
|---|---:|---|
| `reproduce.py` | 16618 | `096de42b1a85f7a27c5fa046a9708d3a1b8b4b98a4304624aabf4e2928d3eb2f` |
| `reproduction/__init__.py` | 72 | `ed9010e9708a159da43e98f7993611d2bf5ef32161c2f3c5381f6ee294ddd199` |
| `reproduction/common.py` | 13098 | `170f6eb17afd40a31538064de15d734c14aa44065978a5f374b5c34811f3444b` |
| `reproduction/figure1.py` | 24311 | `08b93d844eaddd5524ae9be49862f7fff4b668839a8a9dde64bf59e8e2119220` |
| `reproduction/figure1_public.py` | 12261 | `30f23301e2d63630a888eca59847768b475c8ca004362e3a2a85afe25303d688` |
| `reproduction/figure2.py` | 45544 | `fb8e518cfcd4de6ecef4b6965f0a9ef772efb0b7ef4e6593a35df74cb7f51d8f` |
| `reproduction/figure3.py` | 72159 | `449c7ccc6ede6cd6c84b1cb1c42ff149786b7678a964e442f626920d53064c6c` |
| `reproduction/figure4.py` | 31939 | `0dddfa0dbdf6b0a94a59197b1848fc10b54b76241377e6631685fadc1db24e22` |
| `reproduction/figure4ef.py` | 48018 | `ec152c36518d284b143e9dd987054b4e7a0b0c929c806299f4c822df61cbfe27` |
| `reproduction/figure4g.py` | 4986 | `a09beffeb4fb4fe64537c5a3149a761c04eb22948f568a61f15d0f36452e3d62` |
| `reproduction/report.py` | 39756 | `300f88e1e49be20c16014c6d0285277f8c996f695627eeac84f41f3eb63f39a4` |

## Panel results

| Panel | Inputs | AlphaGenome regime | Scored units | API requests | Comparison |
|---|---|---|---:|---:|---|
| 4G | AlphaGenome API (single variant, RNA_SEQ scorer) | ALL_FOLDS | 1 | 1 | PASS |

## Timings

| Step | Status | Seconds |
|---|---|---:|
| 4G: score rs12740374 RNA_SEQ across tissue ontologies | complete | 1.522 |
| 4G: aggregate per-context mean RNA delta and derived adipose row | complete | 0.015 |

## Downloads

No external datasets were downloaded for the selected panels.

## Generated-file manifest

(Omitted here: 7 lines listing every generated file with its SHA-256 -- see run.json for the full audit trail. Available in the actual run directory this report was generated from.)

## Interpretation

A PASS means identities and counts were exact and numerical outputs met the explicit panel-specific equivalence thresholds recorded in `audit/comparison.json`. Figure 2F and non-live comparisons retain `rtol=1e-5, atol=1e-6`; live sequence/ISM panels additionally allow bounded sub-panel-unit API drift while requiring very high reference correlation. The observed maximum differences and plotted-correlation changes remain fully reported. The frozen tables were not used to generate any result.
