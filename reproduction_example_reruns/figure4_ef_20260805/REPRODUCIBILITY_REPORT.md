# SORT1 clean-room reproducibility report

- Run status: **complete**
- Reference comparison: **PASS**
- Started (UTC): `2026-08-05T16:03:16.526895+00:00`
- Repository commit: `d70c3a8e0f839b37eba61defa10cb57eb3713075`
- Repository dirty during run: `True`
- Run directory: `/Volumes/T7/alphaGenome/repro_crash_test/runs/figure4_ef_20260805`
- AlphaGenome credential available during run: `True` (value never recorded)
- Timed step total: `12761.9 s`
- Timed step total across retained attempts: `12761.9 s`
- Peak RSS reported by process: `3911417856 bytes`
- Run-directory size: `4622622695 bytes`

## Scope

The analysis stage is isolated from `outputs/source_data/`. Published source tables are opened only after a successful analysis, by the comparison command.

## Post-run comparison code

- Path: `reproduction/report.py`
- Bytes: `38810`
- SHA-256: `d19b324cb9641c7296a341e5b266a97e75263cb01c3d51545a649392836290e8`

## Executed reproduction code

| Path | Bytes | SHA-256 |
|---|---:|---|
| `reproduce.py` | 16493 | `4c6388a2bc1c06af81b0ffe10f67d64831aebd0b4177cd437de3c660f644ff20` |
| `reproduction/__init__.py` | 72 | `ed9010e9708a159da43e98f7993611d2bf5ef32161c2f3c5381f6ee294ddd199` |
| `reproduction/common.py` | 13098 | `170f6eb17afd40a31538064de15d734c14aa44065978a5f374b5c34811f3444b` |
| `reproduction/figure1.py` | 24311 | `08b93d844eaddd5524ae9be49862f7fff4b668839a8a9dde64bf59e8e2119220` |
| `reproduction/figure1_public.py` | 12261 | `30f23301e2d63630a888eca59847768b475c8ca004362e3a2a85afe25303d688` |
| `reproduction/figure2.py` | 45544 | `fb8e518cfcd4de6ecef4b6965f0a9ef772efb0b7ef4e6593a35df74cb7f51d8f` |
| `reproduction/figure3.py` | 72159 | `449c7ccc6ede6cd6c84b1cb1c42ff149786b7678a964e442f626920d53064c6c` |
| `reproduction/figure4.py` | 31939 | `0dddfa0dbdf6b0a94a59197b1848fc10b54b76241377e6631685fadc1db24e22` |
| `reproduction/figure4ef.py` | 48018 | `ec152c36518d284b143e9dd987054b4e7a0b0c929c806299f4c822df61cbfe27` |
| `reproduction/report.py` | 36808 | `f31c5d920425c5b3422d40d1128fe086a0ebb1eb5b087ad0ebc0e0db1e6105a0` |

## Panel results

| Panel | Inputs | AlphaGenome regime | Scored units | API requests | Comparison |
|---|---|---|---:|---:|---|
| 4E | GRCh38 + GENCODE v46 + 4DN HepG2 Hi-C + AlphaGenome API | ALL_FOLDS | 13517 | 106 | PASS |
| 4F | GRCh38 + GENCODE v46 + 4DN HepG2 Hi-C + AlphaGenome API (shares 4E's run) | ALL_FOLDS | 0 | 0 | PASS |

## Timings

| Step | Status | Seconds |
|---|---|---:|
| 4E/4F: chr1 Hi-C-guided distal-contact 315bp T/G transfer benchmark | complete | 6681.474 |
| 4B/4C: download UCSC-sourced hg38 (the build the working archive actually used) | complete | 117.991 |
| 4C: download GENCODE v46 annotation feather | complete | 10.691 |
| 4E/4F: build chr1 GENCODE promoter/exon catalogue | complete | 1.296 |
| 4E/4F: query observed HepG2 Hi-C contacts for 2066 chr1 protein-coding promoters | complete | 2129.212 |
| 4E/4F: annotate sequence context and select high/low contact-site triples | complete | 90.359 |
| 4B/4C: download UCSC-sourced hg38 (the build the working archive actually used) | complete | 4.238 |
| 4E/4F: build 7-state-per-promoter native/T/G transfer design | complete | 22.128 |
| 4E/4F: score states 1-128 of 13517 | complete | 32.577 |
| 4E/4F: score states 129-256 of 13517 | complete | 35.390 |
| 4E/4F: score states 257-384 of 13517 | complete | 33.425 |
| 4E/4F: score states 385-512 of 13517 | complete | 32.009 |
| 4E/4F: score states 513-640 of 13517 | complete | 31.855 |
| 4E/4F: score states 641-768 of 13517 | complete | 31.234 |
| 4E/4F: score states 769-896 of 13517 | complete | 31.904 |
| 4E/4F: score states 897-1024 of 13517 | complete | 33.538 |
| 4E/4F: score states 1025-1152 of 13517 | complete | 31.085 |
| 4E/4F: score states 1153-1280 of 13517 | complete | 31.479 |
| 4E/4F: score states 1281-1408 of 13517 | complete | 32.477 |
| 4E/4F: score states 1409-1536 of 13517 | complete | 32.081 |
| 4E/4F: score states 1537-1664 of 13517 | complete | 31.006 |
| 4E/4F: score states 1665-1792 of 13517 | complete | 32.979 |
| 4E/4F: score states 1793-1920 of 13517 | complete | 32.921 |
| 4E/4F: score states 1921-2048 of 13517 | complete | 32.407 |
| 4E/4F: score states 2049-2176 of 13517 | complete | 32.372 |
| 4E/4F: score states 2177-2304 of 13517 | complete | 33.383 |
| 4E/4F: score states 2305-2432 of 13517 | complete | 32.726 |
| 4E/4F: score states 2433-2560 of 13517 | complete | 32.622 |
| 4E/4F: score states 2561-2688 of 13517 | complete | 31.986 |
| 4E/4F: score states 2689-2816 of 13517 | complete | 32.691 |
| 4E/4F: score states 2817-2944 of 13517 | complete | 32.552 |
| 4E/4F: score states 2945-3072 of 13517 | complete | 33.089 |
| 4E/4F: score states 3073-3200 of 13517 | complete | 31.857 |
| 4E/4F: score states 3201-3328 of 13517 | complete | 31.435 |
| 4E/4F: score states 3329-3456 of 13517 | complete | 31.369 |
| 4E/4F: score states 3457-3584 of 13517 | complete | 32.161 |
| 4E/4F: score states 3585-3712 of 13517 | complete | 32.367 |
| 4E/4F: score states 3713-3840 of 13517 | complete | 33.618 |
| 4E/4F: score states 3841-3968 of 13517 | complete | 33.046 |
| 4E/4F: score states 3969-4096 of 13517 | complete | 32.774 |
| 4E/4F: score states 4097-4224 of 13517 | complete | 33.181 |
| 4E/4F: score states 4225-4352 of 13517 | complete | 32.296 |
| 4E/4F: score states 4353-4480 of 13517 | complete | 32.297 |
| 4E/4F: score states 4481-4608 of 13517 | complete | 32.901 |
| 4E/4F: score states 4609-4736 of 13517 | complete | 33.491 |
| 4E/4F: score states 4737-4864 of 13517 | complete | 33.187 |
| 4E/4F: score states 4865-4992 of 13517 | complete | 32.586 |
| 4E/4F: score states 4993-5120 of 13517 | complete | 32.765 |
| 4E/4F: score states 5121-5248 of 13517 | complete | 32.927 |
| 4E/4F: score states 5249-5376 of 13517 | complete | 33.806 |
| 4E/4F: score states 5377-5504 of 13517 | complete | 33.459 |
| 4E/4F: score states 5505-5632 of 13517 | complete | 33.386 |
| 4E/4F: score states 5633-5760 of 13517 | complete | 33.051 |
| 4E/4F: score states 5761-5888 of 13517 | complete | 33.800 |
| 4E/4F: score states 5889-6016 of 13517 | complete | 34.254 |
| 4E/4F: score states 6017-6144 of 13517 | complete | 33.400 |
| 4E/4F: score states 6145-6272 of 13517 | complete | 33.304 |
| 4E/4F: score states 6273-6400 of 13517 | complete | 34.122 |
| 4E/4F: score states 6401-6528 of 13517 | complete | 37.064 |
| 4E/4F: score states 6529-6656 of 13517 | complete | 36.676 |
| 4E/4F: score states 6657-6784 of 13517 | complete | 34.085 |
| 4E/4F: score states 6785-6912 of 13517 | complete | 33.907 |
| 4E/4F: score states 6913-7040 of 13517 | complete | 34.203 |
| 4E/4F: score states 7041-7168 of 13517 | complete | 33.868 |
| 4E/4F: score states 7169-7296 of 13517 | complete | 34.672 |
| 4E/4F: score states 7297-7424 of 13517 | complete | 34.501 |
| 4E/4F: score states 7425-7552 of 13517 | complete | 36.652 |
| 4E/4F: score states 7553-7680 of 13517 | complete | 34.640 |
| 4E/4F: score states 7681-7808 of 13517 | complete | 37.220 |
| 4E/4F: score states 7809-7936 of 13517 | complete | 38.198 |
| 4E/4F: score states 7937-8064 of 13517 | complete | 34.443 |
| 4E/4F: score states 8065-8192 of 13517 | complete | 34.519 |
| 4E/4F: score states 8193-8320 of 13517 | complete | 34.983 |
| 4E/4F: score states 8321-8448 of 13517 | complete | 35.295 |
| 4E/4F: score states 8449-8576 of 13517 | complete | 35.292 |
| 4E/4F: score states 8577-8704 of 13517 | complete | 36.163 |
| 4E/4F: score states 8705-8832 of 13517 | complete | 36.842 |
| 4E/4F: score states 8833-8960 of 13517 | complete | 35.532 |
| 4E/4F: score states 8961-9088 of 13517 | complete | 35.484 |
| 4E/4F: score states 9089-9216 of 13517 | complete | 35.802 |
| 4E/4F: score states 9217-9344 of 13517 | complete | 35.918 |
| 4E/4F: score states 9345-9472 of 13517 | complete | 35.837 |
| 4E/4F: score states 9473-9600 of 13517 | complete | 35.867 |
| 4E/4F: score states 9601-9728 of 13517 | complete | 35.668 |
| 4E/4F: score states 9729-9856 of 13517 | complete | 35.859 |
| 4E/4F: score states 9857-9984 of 13517 | complete | 37.018 |
| 4E/4F: score states 9985-10112 of 13517 | complete | 39.664 |
| 4E/4F: score states 10113-10240 of 13517 | complete | 36.136 |
| 4E/4F: score states 10241-10368 of 13517 | complete | 38.403 |
| 4E/4F: score states 10369-10496 of 13517 | complete | 36.741 |
| 4E/4F: score states 10497-10624 of 13517 | complete | 37.312 |
| 4E/4F: score states 10625-10752 of 13517 | complete | 36.204 |
| 4E/4F: score states 10753-10880 of 13517 | complete | 35.943 |
| 4E/4F: score states 10881-11008 of 13517 | complete | 38.630 |
| 4E/4F: score states 11009-11136 of 13517 | complete | 36.342 |
| 4E/4F: score states 11137-11264 of 13517 | complete | 37.720 |
| 4E/4F: score states 11265-11392 of 13517 | complete | 37.253 |
| 4E/4F: score states 11393-11520 of 13517 | complete | 36.389 |
| 4E/4F: score states 11521-11648 of 13517 | complete | 37.137 |
| 4E/4F: score states 11649-11776 of 13517 | complete | 63.192 |
| 4E/4F: score states 11777-11904 of 13517 | complete | 36.634 |
| 4E/4F: score states 11905-12032 of 13517 | complete | 38.966 |
| 4E/4F: score states 12033-12160 of 13517 | complete | 36.964 |
| 4E/4F: score states 12161-12288 of 13517 | complete | 37.317 |
| 4E/4F: score states 12289-12416 of 13517 | complete | 37.816 |
| 4E/4F: score states 12417-12544 of 13517 | complete | 38.091 |
| 4E/4F: score states 12545-12672 of 13517 | complete | 38.612 |
| 4E/4F: score states 12673-12800 of 13517 | complete | 37.777 |
| 4E/4F: score states 12801-12928 of 13517 | complete | 37.181 |
| 4E/4F: score states 12929-13056 of 13517 | complete | 38.045 |
| 4E/4F: score states 13057-13184 of 13517 | complete | 40.075 |
| 4E/4F: score states 13185-13312 of 13517 | complete | 38.494 |
| 4E/4F: score states 13313-13440 of 13517 | complete | 38.280 |
| 4E/4F: score states 13441-13517 of 13517 | complete | 27.118 |
| 4E/4F: compute promoter-level T-G high-vs-low interaction and bootstrap summaries | complete | 3.276 |

## Downloads

| URL | Bytes | SHA-256 | Reused |
|---|---:|---|---|
| https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz | 983659424 | `c1dd87068c254eb53d944f71e51d1311964fce8de24d6fc0effc9c61c01527d4` | False |
| https://storage.googleapis.com/alphagenome/reference/gencode/hg38/gencode.v46.annotation.gtf.gz.feather | 333040258 | `7b10f643d96e1142ef058d9c08487f4360552cea29c8339f195a5a96489dbb4c` | False |

## Generated-file manifest

(Omitted here: 13,581 lines listing every generated file with its SHA-256 -- see run.json for the full audit trail. Available in the actual run directory this report was generated from.)

## Interpretation

A PASS means identities and counts were exact and numerical outputs met the explicit panel-specific equivalence thresholds recorded in `audit/comparison.json`. Figure 2F and non-live comparisons retain `rtol=1e-5, atol=1e-6`; live sequence/ISM panels additionally allow bounded sub-panel-unit API drift while requiring very high reference correlation. The observed maximum differences and plotted-correlation changes remain fully reported. The frozen tables were not used to generate any result.
