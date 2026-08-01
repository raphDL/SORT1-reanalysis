# Analysis code

Analysis scripts generate compact source tables from public inputs and/or
AlphaGenome predictions. They may be computationally expensive and may require
API credentials. They must not write directly into the working manuscript
folders.

Each analysis entry point must:

1. read paths and constants from `config.yaml` or explicit command-line
   arguments;
2. record the model regime, track identifiers, genome build, code revisions,
   and input checksums;
3. write a compact panel source table into `outputs/source_data/`;
4. write a run manifest into `outputs/run_manifests/`;
5. never overwrite a frozen output without retaining the prior version.
