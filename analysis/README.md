# Analysis code

Analysis scripts generate compact source tables from public inputs and/or
AlphaGenome predictions. They may be computationally expensive and may require
API credentials. Public clean-room reproduction uses `reproduce.py` and writes
only to an isolated run directory. Frozen `outputs/source_data/` files are
reference answers, not inputs to the analysis stage.

Each analysis entry point must:

1. read paths and constants from `config.yaml` or explicit command-line
   arguments;
2. record the model regime, track identifiers, genome build, code revisions,
   and input checksums;
3. write new predictions and compact panel tables only inside the selected run
   directory;
4. write timing, environment, input-checksum, API-call, and output-checksum
   records under that run directory's `audit/` folder;
5. leave `outputs/source_data/` untouched until the separate comparison stage,
   which reads it without modifying it.
