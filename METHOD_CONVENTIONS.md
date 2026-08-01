# Coordinate, allele and score conventions

This file records cross-analysis conventions that otherwise become difficult
to recover from individual plotting scripts.

## Reference genome and rs12740374

- AlphaGenome sequence analyses use GRCh38.
- rs12740374 is represented as `chr1:109274968 G>T` on the forward strand.
- `REF` and `ALT` denote the GRCh38 reference G allele and alternate/minor T
  allele. They are not labeled wild type and mutant.
- Unless stated otherwise, an allele-specific delta is `ALT - REF` (T minus
  G). Population effect estimates are explicitly reoriented before comparison
  rather than assumed to share this axis.

## Population tagging analyses

- Tagging coefficients were derived from phased 1000 Genomes Phase 3 European
  haplotypes.
- For variant *i*, the allele-oriented coefficient is
  `t_i = Cov(G_i, G_rs12740374) / Var(G_i)` after aligning both genotype
  vectors to the reported allele axis.
- The explanatory tagging model is
  `M_i = D_i + t_i D_rs12740374`; rs12740374 itself receives its direct score
  only once. This construction assumes rs12740374 is the driver and is not an
  independent causal test.

## Wang repair-product coordinates

- Wang et al. repair coordinates were reported on hg18 and mapped to the local
  GRCh38 sequence with a verified offset of -344,145 bp. The aligned local
  reference matched 379 of 380 positions.
- The inferred Cas9 cleavage boundary is
  `chr1:109274966|109274967` (GRCh38).
- Insertions and deletions were reconstructed on the rs12740374-T background.
  The model input length was held fixed by trimming or extending only at the
  distal downstream boundary.

## RNA and chromatin endpoints

- Gene-level RNA summaries use strand-compatible tracks and windows centered
  on the annotated transcription start site. The exact model context and
  output window for each panel are recorded in
  `outputs/run_manifests/panel_analysis_parameters.tsv`.
- Figure 1C uses the recommended exon-masked RNA variant scorer and a
  1,048,576-bp input. Figure 2 repair analyses use the postnatal liver RNA
  target designated track 4 in the stored output metadata and a 524,288-bp
  input.
- Kircher regulatory predictions use 16,384-bp windows; Kircher RNA
  predictions use 131,072-bp windows. These are not interchangeable with the
  locus-scale sequence contexts.

## Contact maps

- Experimental HepG2 Hi-C values are observed/expected (O/E) contacts, so
  genomic-distance decay has already been normalized by the source matrix.
- Virtual-4C profiles are anchor rows/slices extracted from a two-dimensional
  contact matrix; they are derived summaries, not a separate experimental or
  AlphaGenome output type.
- Experimental and predicted contact values undergo different preprocessing.
  Their absolute scales are not treated as calibrated; comparisons use
  percentiles among contacts spanning the same genomic distance.

## Scramble and transfer designs

- Scramble retention is based on matched G and T constructs and normalized to
  the intact-locus G-to-T effect.
- All module-transfer and distal-contact experiments use length-preserving
  replacement, not sequence insertion: the 315-bp donor replaces 315 bp of the
  recipient sequence so model input coordinates and total length remain fixed.
