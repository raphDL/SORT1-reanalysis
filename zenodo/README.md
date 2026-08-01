# Zenodo data deposit

The public data deposit consists of the four large derived AlphaGenome tables
listed in
`outputs/run_manifests/zenodo_pending_large_outputs.tsv` (736,108,954 bytes in
total). Each file has a SHA-256 digest and an exact regeneration-script
pointer. The files are deliberately not duplicated in git.

To stage a verified upload directory from the frozen working archive:

```bash
python analysis/stage_zenodo_deposit.py \
  --workspace-root /path/to/alphaGenome \
  --output-dir /path/to/zenodo-upload
```

After creating the Zenodo record, update
`outputs/run_manifests/zenodo_deposit_status.tsv`, `README.md`, and
`CITATION.cff` with the DOI and record URL. The strict release gate rejects a
pending status or DOI.

The four files are derived model outputs, not redistributed public input
datasets. Confirm AlphaGenome prediction-redistribution terms before making
the Zenodo record public.
