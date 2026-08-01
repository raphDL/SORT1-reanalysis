# Prespecified plan: held-out replication of the rs12740374 tissue profile

Date frozen: 2026-07-26

## Question

Does the tissue pattern shown in Figure 4A with the default `ALL_FOLDS`
AlphaGenome model reproduce when the locus is scored with `FOLD_0`, for which
the rs12740374 region is held out from training?

## Frozen design

- Variant: hg38 `chr1:109274968:G>T` (rs12740374).
- Sequence window: 1,048,576 bp, identical to the source analysis.
- Readout: the recommended RNA-seq variant scorer.
- Genes: SORT1, PSRC1, and CELSR2.
- Tissue groups: liver, adipose (mean of subcutaneous and visceral adipose),
  spleen, lung, brain, heart, and kidney.
- Track aggregation: reproduce the source script exactly within each ontology
  and gene, then average the two adipose ontology estimates.
- Comparator: the frozen `ALL_FOLDS` source table already used for Figure 4A.

## Frozen evaluation

Primary descriptive checks:

1. The FOLD_0 liver estimate is positive for each of the three genes.
2. Liver has the largest positive estimate among the seven displayed tissue
   groups for each gene.

Secondary stability checks:

- Pearson and Spearman correlation across all 21 tissue-by-gene cells.
- Sign agreement across the 21 cells.
- Rank of liver within the seven tissues for each gene.
- A nonparametric bootstrap 95% confidence interval for Spearman correlation,
  resampling the 21 tissue-by-gene cells with seed 12740374.

No post hoc threshold will convert this analysis into a pass/fail claim.
The manuscript language will report the observed estimates and uncertainty.

## Placement rule

- If both primary checks hold, the FOLD_0 result can support the Figure 4A
  tissue-profile statement in the main text.
- If either primary check fails, Figure 4A remains an `ALL_FOLDS` exploratory
  analysis and the held-out result is reported as a limitation or supplementary
  result.

