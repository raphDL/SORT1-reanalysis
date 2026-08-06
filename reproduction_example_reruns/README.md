# Example clean-room re-runs

This directory holds the audit trail from real executions of
`reproduce.py` (see the top-level [README](../README.md) and
[`reproduce.py`](../reproduce.py) / [`reproduction/`](../reproduction/)) —
evidence that the reproducibility capsule was actually run end-to-end, not
just written. Each subdirectory is `run.json` (full machine-readable audit:
every step, timing, download, and API call count) plus
`REPRODUCIBILITY_REPORT.md` (the human-readable summary `reproduce.py`
itself generates), copied from a completed run directory. The very large
per-generated-file checksum manifest section of each report is trimmed for
size; `run.json` still contains the full record.

These runs were executed on a separate machine/volume from the
development checkout, each starting from empty `raw/`, `predictions/`, and
`derived/` directories — so every download and every AlphaGenome
prediction in them was freshly fetched/scored for that run, not reused
from anywhere else. See `REPRODUCIBILITY_NEXT_STEPS.md` at the repository
root for the full narrative (bugs found and fixed, verification steps,
what's not yet done).

| Run | Panels | Started (UTC) | Result | Repository commit |
|---|---|---|---|---|
| `figure1_full_public_20260802T224000Z` | 1B, 1C, 1D, 1E, 1F | 2026-08-02T20:35:59Z | PASS | `f0423a7` |
| `figure2_public_inputs_20260802T192638Z` | 2B, 2C, 2E, 2F | 2026-08-02T19:56:48Z | PASS | `f0423a7` |
| `figure3_all_panels_20260805` | 3A, 3B, 3C, 3E, 3F, 3G | 2026-08-05T11:26:31Z | PASS | `23171f0` |
| `figure4_bc_20260805` | 4B, 4C | 2026-08-05T12:40:48Z | PASS | `81bb2bf` |
| `figure4_ef_20260805` | 4E, 4F | 2026-08-05T16:03:16Z | PASS | `d70c3a8` |
| `figure4g_20260805` | 4G | 2026-08-05T18:24:34Z | PASS | `328cd24` |
| `figure4h_20260805_stopped_early` | 4H | 2026-08-05T21:16:44Z | stopped early (256/~900,000 predictions; not a failure — see note below) | `99e9931` |

Panels 4A, 4D, and 4I are non-computational author-layout schematics and
are out of scope for `reproduce.py`. Panel 4J was scoped (see
`REPRODUCIBILITY_NEXT_STEPS.md`) but its run was never started.

The "Repository commit" column is the exact commit each run executed
against — later commits in this repository's history than that hash are
documentation/wiring for panels done in subsequent runs, not changes that
would alter an already-passing panel's result.

`figure4h_20260805_stopped_early` was deliberately halted mid-run by
project decision (the wall-clock/API cost for the remaining ~899,700
predictions wasn't judged worth spending, especially after a targeted
smoke-test finding — documented in `REPRODUCIBILITY_NEXT_STEPS.md` —
that individual-SNP effects at this panel's resolution sit at or below
AlphaGenome's own measurement drift over time). It is included here for
transparency about exactly how far that run got, not as a passing result.
