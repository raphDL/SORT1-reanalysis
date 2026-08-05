# Data

`SOURCES.tsv` records every public or restricted input used by the selected
analyses, including the exact release or accession, the release-audit access
date, redistribution status, and a checksum whenever the original artifact
was retained. `reproduce.py` automates public downloads for its supported
clean-room panels and records their observed hashes in the run directory.
Other inputs are very large, access-controlled, or governed by
source-specific terms; for those, the stable accession/URL remains the
retrieval recipe. Frozen compact derivatives under `outputs/source_data/` are
publication references for post-run comparison, not substitutes for raw
inputs during clean-room analysis.

No public dataset should be committed merely for convenience. If
redistribution rights are uncertain, include accession, checksum, and a fetch
or preparation recipe instead.
