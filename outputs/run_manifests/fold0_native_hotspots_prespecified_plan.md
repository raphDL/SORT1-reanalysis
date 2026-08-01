# Prespecified plan: held-out audit of native 100-kb SNV hotspots

Date frozen: 2026-07-26

## Question

Do the strongest coordinated native SNV effects in the existing 100-kb
`ALL_FOLDS` liver, memory-T-cell, and CD14-positive-monocyte scans retain their
direction under `FOLD_0`?

## Frozen selection

- Source interval: hg38 chr1:109,224,968-109,324,968.
- Source statistic: `synergy_score`, defined as the minimum of the three
  positive RNA deltas when SORT1, PSRC1, and CELSR2 are all increased, and zero
  otherwise.
- Select the top five variants by `synergy_score` separately in liver,
  CD8-positive memory T cell, and CD14-positive monocyte.
- Selection is made exclusively from the existing `ALL_FOLDS` tables before
  any FOLD_0 score is generated.
- Exact REF and ALT alleles are retained. Duplicate variants across contexts
  remain separate context-specific tests.

## Frozen scoring

- Model: `FOLD_0`.
- Prediction window: 524,288 bp centered on rs12740374.
- Readout: RNA-seq signal summed over +/-2 kb of the annotated TSS for SORT1,
  PSRC1, and CELSR2.
- Tracks:
  - liver: `UBERON:0002107`, polyA-plus RNA-seq;
  - T cell: `CL:0000909`, polyA-plus RNA-seq;
  - monocyte: `CL:0001054`, polyA-plus RNA-seq.
- For each context, match the same metadata-defined track used by the source
  scan.

## Frozen evaluation

- Primary: proportion of the 15 context-specific candidates that remain
  positive for all three genes under FOLD_0.
- Secondary: sign agreement for all 45 gene-level effects; Pearson and
  Spearman correlation of gene-level effects; Spearman correlation of the five
  within-context `synergy_score` ranks.
- Exact binomial 95% confidence intervals will accompany sign-agreement
  proportions.

## Interpretation

This is a selected-candidate replication audit, not an unbiased estimate of
genome-wide performance. It can support the direction of named hotspots but
cannot validate the dense atlas or its complete rank ordering.

