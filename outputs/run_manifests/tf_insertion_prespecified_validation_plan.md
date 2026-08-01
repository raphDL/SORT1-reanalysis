# Prespecified validation plan for the Figure 4 motif-replacement scan

Date frozen: 2026-07-26

This plan was written before generating or inspecting matched-shuffle or
AlphaGenome FOLD_0 results. Its purpose is to prevent decisions about main-text
placement from depending on how the validation results land.

## Frozen discovery data

The candidate sites will be selected only from:

`dense_100kb_8tf_3tracks_variant_summary.tsv`

SHA-256:

`16d6e6b0d2f056b93ef58d67fcdf24971efb1887cec6049a0a32da1c885eb941`

The corresponding design table has SHA-256:

`57dfce1ee5d11d9738cb58a4d61a7dbc3daf7f413090642f9472bd97be45238d`

No matched-shuffle or FOLD_0 result will be used to select candidate sites.

## Corrected discovery score

1. Use raw AlphaGenome RNA deltas; do not subtract the rs12740374 effect.
2. At each motif and position, retain the orientation with the largest
   **positive** effect, not the largest absolute effect.
3. Treat the sum across SORT1, PSRC1, and CELSR2 as **net locus RNA output**.
4. Call a result **coordinated activation** only when all three gene-specific
   deltas are positive. The coordinated score is then the smallest of the
   three positive deltas; otherwise it is zero.

## Frozen primary candidates

For each cellular context, select the five highest raw net-output sites among
the prespecified context-associated motifs:

- liver: HNF1B, HNF4A, FOXA2, and CEBP;
- T cell: SPI1, RUNX1, IRF8, and RELA;
- CD14-positive monocyte: SPI1, RUNX1, IRF8, and RELA.

Selected sites must be separated by at least 500 bp. The highest-scoring site
is retained first, followed iteratively by the next site at least 500 bp from
all retained sites. This produces 15 primary candidates.

For secondary analysis, retain the three highest sites for every
motif-by-context combination using the same 500-bp separation rule. Secondary
results cannot rescue failure of the primary analysis.

## Matched sequence controls

For every frozen candidate:

1. Generate up to 199 unique mononucleotide-composition-preserving shuffles of
   the corresponding motif sequence using a fixed recorded seed.
2. Replace the identical native interval, so motif and control constructs have
   the same length and remove the same reference bases.
3. Score both orientations for every shuffled sequence and retain the largest
   positive orientation using the same rule applied to the motif.
4. Define the control-adjusted effect as the motif's raw effect minus the
   median matched-shuffle effect.
5. Calculate a one-sided empirical P value as
   `(1 + number of shuffle effects >= motif effect) / (1 + number of shuffles)`.
6. Apply Benjamini-Hochberg correction across the 15 primary candidates.

If fewer than 199 unique shuffles exist for a motif, enumerate all unique
non-native, non-reverse-complement arrangements and report the resulting
empirical resolution.

## Held-out-locus replication

The 15 frozen primary motif constructs and all of their matched controls will
be rescored with AlphaGenome FOLD_0. Candidate selection and thresholds will
not be changed after inspecting these results.

## Annotation concordance

Primary candidates will be compared with context-matched chromatin
accessibility and H3K27ac annotations. Enrichment will be evaluated against
position-matched null sites that preserve distance to the nearest analyzed TSS
and baseline accessibility strata. Annotation is supporting evidence and will
not be described as validation of transcription-factor binding.

## Main-text stopping rule

The motif-replacement panel remains a principal main-text result only if all
of the following are met:

1. **Matched controls:** at least 8 of 15 primary candidates have a positive
   control-adjusted effect and empirical FDR `q < 0.10`, with at least one
   passing candidate in each cellular context.
2. **FOLD_0:** at least 10 of 15 retain a positive control-adjusted effect;
   the all-fold versus FOLD_0 control-adjusted effects have Spearman
   correlation `rho >= 0.50`; and the top candidate in each context retains a
   positive effect with empirical `P < 0.05`.
3. **Observed annotations:** the primary candidate set shows at least
   twofold enrichment in matched accessibility or H3K27ac annotations with
   permutation `P < 0.05`, and each context has at least one supported
   candidate.
4. **Claim-specific requirement:** any candidate described as coordinated
   activation has positive gene-specific effects for SORT1, PSRC1, and CELSR2
   in FOLD_0.

If criterion 1 or 2 fails, the motif-replacement analysis will be demoted to
the supplement and described only as exploratory model-based design. If only
criterion 3 fails, it will also be demoted unless the main-text claim is
restricted explicitly to synthetic sequence receptivity rather than native
regulatory hotspots. Failure of criterion 4 removes the coordinated-activation
claim but does not invalidate a clearly labelled net-output result.

## Prohibited post hoc changes

After matched-shuffle results are inspected, do not:

- change candidate sites or their separation distance;
- switch between raw sum, minimum-positive, maximum-gene, or
  rs12740374-normalized scores to obtain a passing result;
- change orientation selection;
- drop an unfavorable motif or cellular context;
- replace the empirical control family;
- redefine the main-text stopping thresholds.

Any additional analysis must be labelled exploratory.
