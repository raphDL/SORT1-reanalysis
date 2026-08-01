# Prespecified plan: genome-wide liver caQTL benchmark

Date frozen: 2026-07-26

## Question

Can AlphaGenome distinguish experimentally significant liver chromatin-
accessibility QTLs from matched tested variants and recover their allelic
direction beyond the SORT1 locus?

## Source data

- Currin et al. liver caQTL lead variants from the 1-kb analysis.
- Chromosome-wide nominal association statistics from the released archive.
- hg38 coordinates and alleles as released by the authors.

## Frozen eligibility

- Autosomal biallelic SNVs only.
- REF/ALT must match the hg38 reference representation in the released variant
  identifier.
- Minor-allele frequency at least 0.05.
- Imputation R2 at least 0.8 when the field is available.
- One row per unique variant in the primary benchmark. When a significant lead
  variant is reported for multiple peaks, retain the variant-peak association
  with the smallest q value (then smallest nominal p value as a deterministic
  tie-break). Controls are also unique by variant; if a tested variant is
  associated with multiple eligible null peaks, retain one matched pair at
  random under the frozen seed. This prevents one AlphaGenome variant score
  from being counted repeatedly.

## Frozen sample

- Positives: 500 significant 1-kb lead caQTL SNVs, sampled with seed 12740374
  after stratification by chromosome.
- Controls: 500 tested SNV-peak pairs with nominal p >= 0.5, sampled without
  replacement from the same chromosomes and matched to positives by:
  - MAF bin: [0.05,0.10), [0.10,0.20), [0.20,0.30), [0.30,0.40), [0.40,0.50];
  - absolute distance to peak center: [0,100), [100,250), [250,500),
    [500,1000] bp;
  - transition versus transversion.
- A control variant cannot be a significant lead caQTL or rs12740374.
- If an exact stratum lacks a control, matching relaxes first by transition
  class and then to an adjacent MAF bin; every relaxation is recorded.

The sampled manifest and its SHA256 digest will be written before any
AlphaGenome result is requested.

## Frozen prediction

- Model: `ALL_FOLDS`.
- Primary output: ATAC-seq `DIFF_LOG2_SUM` in liver/HepG2-associated tracks,
  using the same aggregation as the existing Currin/SORT1 analysis.
- Variant effect is REF-to-ALT.
- Experimental beta will be reoriented to REF-to-ALT using the released
  effect-allele (`EA`) and non-effect-allele (`NEA`) fields. Rows with
  irreconcilable alleles are excluded before scoring.

## Frozen primary endpoints

1. Classification of positive versus matched control by absolute predicted
   ATAC effect: ROC AUC and precision-recall AUC, with 2,000 paired bootstrap
   replicates (seed 12740374).
2. Direction agreement among positive caQTLs: proportion for which the sign of
   the AlphaGenome REF-to-ALT score matches the sign of the REF-to-ALT
   experimental beta, with exact binomial 95% confidence interval.

Secondary endpoints:

- Spearman correlation between experimental beta and predicted ATAC effect
  among positives.
- Performance stratified by motif-disruption annotation, effect magnitude,
  MAF, and distance to peak center.
- Sensitivity to liver-only versus HepG2-only track aggregation.

## Interpretation

This is an external association benchmark, not a prospective assay and not
proof of causality for every lead variant. LD can make the reported lead a
proxy rather than the molecular driver; classification is therefore expected
to be more robust than effect-size correlation.
