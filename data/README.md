# Data

`SOURCES.tsv` records every public or restricted input used by the selected
analyses, including the exact release or accession, the release-audit access
date, redistribution status, and a checksum whenever the original artifact
was retained. Raw bulk downloads are deliberately not automated here: several
inputs are very large, access-controlled, or governed by source-specific
terms. The stable accession/URL is therefore the retrieval recipe, while the
checksummed compact derivatives used by each panel live in
`outputs/source_data/`.

No public dataset should be committed merely for convenience. If
redistribution rights are uncertain, include accession, checksum, and a fetch
or preparation recipe instead.
