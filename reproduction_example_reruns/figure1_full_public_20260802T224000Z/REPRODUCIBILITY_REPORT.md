# SORT1 clean-room reproducibility report

- Run status: **complete**
- Reference comparison: **PASS**
- Started (UTC): `2026-08-02T20:35:59.409531+00:00`
- Repository commit: `f0423a73ae04fcd50b1cf915ab0d8fa5791e709f`
- Repository dirty during run: `True`
- Run directory: `/Volumes/T7/alphaGenome/repro_crash_test/runs/figure1_full_public_20260802T224000Z`
- AlphaGenome credential available during run: `True` (value never recorded)
- Timed step total: `353.8 s`
- Timed step total across retained attempts: `353.8 s`
- Peak RSS reported by process: `981516288 bytes`
- Run-directory size: `7984594057 bytes`

## Scope

The analysis stage is isolated from `outputs/source_data/`. Published source tables are opened only after a successful analysis, by the comparison command.

## Post-run comparison code

- Path: `reproduction/report.py`
- Bytes: `22078`
- SHA-256: `ec539ebe1012f01d521683becb3647374b4dc0ba5850b9f18cb7715dc95bb676`

## Executed reproduction code

| Path | Bytes | SHA-256 |
|---|---:|---|
| `reproduce.py` | 12822 | `0318fb5c0be41275a147893662c72753805b7d46736e5254764a4aa57ce08f5b` |
| `reproduction/__init__.py` | 72 | `ed9010e9708a159da43e98f7993611d2bf5ef32161c2f3c5381f6ee294ddd199` |
| `reproduction/common.py` | 13098 | `170f6eb17afd40a31538064de15d734c14aa44065978a5f374b5c34811f3444b` |
| `reproduction/figure1.py` | 24311 | `08b93d844eaddd5524ae9be49862f7fff4b668839a8a9dde64bf59e8e2119220` |
| `reproduction/figure1_public.py` | 12157 | `a9626d41daf66048966a1720be67e7dd6700e1621360b6f68d2835ad8d29d195` |
| `reproduction/figure2.py` | 43967 | `9b7bfd442e0bbd6d8cd0ff74e76f49704033281abb22db4d4a4197bf2f8eccde` |
| `reproduction/report.py` | 22017 | `9f07a68be11bccb5a381435c6f9f6b71d7b3a74c3f7c2395e79b7395e2fd3ea3` |

## Panel results

| Panel | Inputs | AlphaGenome regime | Scored units | API requests | Comparison |
|---|---|---|---:|---:|---|
| 1B | AlphaGenome API | ALL_FOLDS | 1 | 1 | PASS |
| 1C | GLGC + GTEx v7 liver + 1000 Genomes EUR + AlphaGenome API | ALL_FOLDS | 111 | 14 | PASS |
| 1D | 4DN HepG2 observed Hi-C | none | 0 | 0 | PASS |
| 1E | AlphaGenome API | FOLD_0 | 1 | 1 | PASS |
| 1F | 4DN HepG2 Hi-C + AlphaGenome contact map | FOLD_0 | 0 | 0 | PASS |

## Timings

| Step | Status | Seconds |
|---|---|---:|
| Figure 1C: stage GTEx v7 and 1000 Genomes public inputs | complete | 9.091 |
| 1B: AlphaGenome locus tracks and rendering | complete | 7.330 |
| 1C-middle: download GLGC, LiftOver chain, and hg38 | complete | 165.912 |
| 1C-middle: reconstruct the 111-variant set | complete | 2.330 |
| 1C-middle: score variants 1-8 | complete | 3.000 |
| 1C-middle: score variants 9-16 | complete | 2.440 |
| 1C-middle: score variants 17-24 | complete | 3.155 |
| 1C-middle: score variants 25-32 | complete | 2.790 |
| 1C-middle: score variants 33-40 | complete | 2.721 |
| 1C-middle: score variants 41-48 | complete | 2.544 |
| 1C-middle: score variants 49-56 | complete | 2.892 |
| 1C-middle: score variants 57-64 | complete | 3.000 |
| 1C-middle: score variants 65-72 | complete | 2.717 |
| 1C-middle: score variants 73-80 | complete | 2.600 |
| 1C-middle: score variants 81-88 | complete | 3.713 |
| 1C-middle: score variants 89-96 | complete | 2.707 |
| 1C-middle: score variants 97-104 | complete | 2.666 |
| 1C-middle: score variants 105-111 | complete | 3.211 |
| 1C: reconstruct GTEx and EUR-tagging panels | complete | 95.283 |
| 1E: AlphaGenome FOLD_0 contact map and rendering | complete | 2.003 |
| 1D/1F: extract public 4DN HepG2 observed Hi-C | complete | 31.684 |

## Downloads

| URL | Bytes | SHA-256 | Reused |
|---|---:|---|---|
| supplied:<local-machine>/alphaGenome/investigation/SORT1_cholesterol_musunuru_2010/Liver.allpairs.txt.gz | 3702606276 | `372afa081939868407afef0322638b8a191e9f05b1451b6a873fdd04e5d81f67` | False |
| supplied:<local-machine>/alphaGenome/investigation/SORT1_cholesterol_musunuru_2010/enhanced_results/haplotype_context/integrated_call_samples_v3.20130502.ALL.panel | 55156 | `b4023dc6ee2d62ee89c8d4d347db4d348e65518d66d346574cdae7a4bbd76858` | False |
| supplied:<local-machine>/alphaGenome/investigation/SORT1_cholesterol_musunuru_2010/enhanced_results/haplotype_context/1kgp_highcov_grch38_sort1_1p13.vcf.gz | 1223305 | `8b5ee6b3c7956cac6af185ce55125f89dd4952da27fc31554e75ca53aeb1f3b0` | False |
| supplied:<local-machine>/alphaGenome/investigation/SORT1_cholesterol_musunuru_2010/enhanced_results/haplotype_context/1kgp_highcov_grch38_sort1_1p13.vcf.gz.tbi | 370 | `e7769ed6b2ddd0c5ba382fc71d89c37918308561f350961c9d4a75b1a5d318f5` | False |
| https://csg.sph.umich.edu/willer/public/lipids2013/jointGwasMc_LDL.txt.gz | 58851380 | `ddf75ed591ffcce8c3e81ae2f4c77f4428cd7dbf4628b24467aa8cf6e85e3dce` | False |
| https://imputationserver.sph.umich.edu/resources/chain/hg19_to_hg38.over.chain.gz | 227698 | `5c0598e500ceb5a78c73086929e8ef993aec309bcafb595139b53d440b125a1d` | False |
| https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz | 872949833 | `fb4243ebb014caf27111f24dd62b7ce42160f28581da6f8fcd6cba5977778d02` | False |
| https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/wfoutput/25104375-a588-46e6-a382-663cee6c332f/4DNFICSTCJQZ.hic | 14732267572 | `S3 ETag ca98fe976d7321969696347d167ba35c-1757` | False |

## Generated-file manifest

(Omitted here: 34 lines listing every generated file with its SHA-256 -- see run.json for the full audit trail. Available in the actual run directory this report was generated from.)

## Interpretation

A PASS means identities and counts were exact and numerical outputs met the explicit panel-specific equivalence thresholds recorded in `audit/comparison.json`. Figure 2F and non-live comparisons retain `rtol=1e-5, atol=1e-6`; live sequence/ISM panels additionally allow bounded sub-panel-unit API drift while requiring very high reference correlation. The observed maximum differences and plotted-correlation changes remain fully reported. The frozen tables were not used to generate any result.
