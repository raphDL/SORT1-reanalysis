# Derived outputs

`source_data/` contains one compact machine-readable table per published
panel. `run_manifests/` contains frozen selections, model assignments,
provenance audits, and SHA-256 manifests. Large prediction tensors and caches
are intentionally excluded from Git and will be deposited separately.

`manuscript_results.tsv` is the machine-generated single source of truth for
headline numerical values quoted in Results, Methods, legends, and the
abstract. It currently contains 89 verified result rows; each has an analysis
identifier and a release source-table pointer. Regenerate it with
`python analysis/build_manuscript_results.py` after changing any source table.

`run_manifests/source_data_provenance.tsv` maps every compact release table to
all byte-identical copies located in the frozen working archive. This is
separate from `source_data_sha256.tsv`, which validates the release copy
itself.
