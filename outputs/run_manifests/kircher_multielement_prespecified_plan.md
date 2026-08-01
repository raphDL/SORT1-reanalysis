# Prespecified plan: multi-element Kircher satMutMPRA benchmark

Date frozen: 2026-07-26

## Feasibility question

Can the published Kircher et al. saturation-mutagenesis measurements for
multiple regulatory elements be mapped unambiguously to hg38 genomic alleles
and to assay-matched AlphaGenome outputs?

The local project currently contains processed measurements only for the SORT1
and SORT1.2 constructs. The complete published source data must therefore be
retrieved and audited before scoring.

## Inclusion criteria

An element is eligible only when all of the following are available:

1. an unambiguous hg38 genomic interval and strand;
2. reference and alternate alleles that can be reconciled with the hg38
   reference sequence;
3. a quantitative experimental effect for individual substitutions;
4. the experimental assay cell type;
5. at least one assay-matched AlphaGenome chromatin output in that context.

Elements failing any criterion will be listed with the reason for exclusion.
No genomic placement will be guessed from construct names.

## Frozen primary comparison

- Unit: matched single-nucleotide substitution.
- Experimental response: published signed reporter effect.
- Predicted response: signed AlphaGenome accessibility effect in the
  assay-matched context, evaluated at the native genomic locus with a 16,384-bp
  prediction interval and `DIFF_LOG2_SUM` center-mask scoring. ATAC is used
  where a matched track exists; DNase is the prespecified accessibility
  fallback.
- Within-element performance: Spearman correlation.
- Across-element summary: random-effects meta-analysis of Fisher-transformed
  correlations, accompanied by the full distribution of element-level
  correlations.

Secondary comparisons may use H3K27ac and RNA only where matched tracks and a
defensible target-gene definition exist. They will not replace accessibility
as the primary endpoint.

## Frozen minimum evidence

- At least five eligible elements, each with at least 100 matched substitutions,
  are required for a “multi-element benchmark.”
- If this minimum is not met, the folder will contain a feasibility report and
  exclusion table, and no generalization claim will be made.

## Context audit and frozen eligible set

The metadata-only AlphaGenome audit identified six unique eligible elements:

| element | experimental cell line | AlphaGenome context | output |
|---|---|---|---|
| F9 | HepG2 | EFO:0001187 | ATAC |
| FOXE1 | HeLa | EFO:0002791 (HeLa-S3) | DNase |
| LDLR | HepG2 | EFO:0001187 | ATAC |
| MYC rs11986220 | LNCaP | EFO:0005726 | DNase |
| PKLR | K562 | EFO:0002067 | ATAC |
| SORT1 | HepG2 | EFO:0001187 | ATAC |

LDLR and LDLR.2 are averaged at matched substitutions. SORT1 and SORT1.2 are
averaged; SORT1-flip is excluded from the primary benchmark because it is an
orientation control, not an independent element. PKLR-48h is used because the
authors' element table specifies the 48-h assay. Only GRCh38 SNVs represented
by at least 10 barcodes are included.

## Interpretation limits

The MPRA tests episomal reporter activity, whereas AlphaGenome predicts native
genomic outputs. Agreement tests shared sequence grammar; disagreement may
reflect chromosomal context rather than model failure.
