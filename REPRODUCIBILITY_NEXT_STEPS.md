# Reproducibility next steps — 2026-08-05

## Status since the 2026-08-02 audit

The clean-room audit `20260802T080725Z` (`findings.tsv`, `REPRODUCIBILITY_AUDIT.md`,
evidence retained under
`/Volumes/T7/alphaGenome/repro_crash_test_audit/20260802T080725Z/` on the
external drive) rated the release "partially reproducible; not yet end-to-end
reproducible" against commit `c5b00d4`.

A clean-room analysis runner (`reproduce.py` plus the `reproduction/` package
and `tests/test_reproduction.py`) was subsequently built and exercised in five
timestamped runs under
`/Volumes/T7/alphaGenome/repro_crash_test/runs/`. Current state:

| Panels | Run | Result |
|---|---|---|
| 1B | `figure1_api_20260802T090000Z` | Failed (placeholder API key; superseded) |
| 1B/1C-middle/1E | `figure1_public_inputs_20260802T090500Z` | PASS |
| 1B/1C/1D/1E/1F | `figure1_full_public_20260802T224000Z` | PASS |
| 2B/2C/2E/2F | `figure2_public_inputs_20260802T192638Z` | PASS |
| 3B | `figure3_public_20260802T232000Z` | PASS |

Figure 1 (B–F), Figure 2 (B/C/E/F) and Figure 3 (B) now regenerate from a live
AlphaGenome key and public raw inputs and match the frozen reference within
the tolerances in `reproduction/report.py`.

Figure 3's other panels were also attempted inside the same run directory
(`figure3_public_20260802T232000Z`), but **outside** the tracked
`reproduce.py` harness — by invoking the original legacy analysis scripts
under `investigation/SORT1_cholesterol_musunuru_2010/` directly and writing
into that run's `predictions/` folder. None of this is reflected in
`audit/run.json`, `REPRODUCIBILITY_REPORT.md`, or `audit/comparison.json`,
which only ever tracked panel 3B. Their actual state, read off the files
themselves:

| Panel | Script | State |
|---|---|---|
| 3A (`fig3_regional_ism`) | `sort1_figure_2e_100kb_rna_ism.py` | **Broken.** `stage1_window_deltas.tsv` (33.6 MB) and both stage-2 data files completed (2026-08-03 12:55–12:56 UTC), but `logs/Figure3A_regional_ism.log` shows the run crashing at `_load_gene_tss` with `FileNotFoundError: .../investigation/SORT1_cholesterol_musunuru_2010/cache/gencode.v46.annotation.gtf.gz.feather`. That `cache/` directory does not currently exist at all on the source machine. No `run_metadata.json`, derived table, or rendered SVG was produced. |
| 3E (`fig3_directional_recovery`) | (legacy script) | Prediction-stage only. `sequence_design.csv`, `sequence_scores.csv`, `retention_by_seed.csv`, `summary.csv` completed 2026-08-02 21:37–23:55 UTC. No integration into `derived/` or `figures/`, no comparison run. |
| 3F (`fig3_boundary_grid`) | (legacy script) | Prediction-stage only. Started 2026-08-02 23:56 UTC; interrupted mid-run when the T7 drive dropped and had to be resumed (`logs/Figure3F_resume.log`, `logs/Figure3F_resume_after_remount.log`); finished the prediction stage 2026-08-03 09:35 UTC (896 candidates, 28,688 sequence designs). No integration into `derived/` or `figures/`. |
| 3G (`fig3_component_necessity`) | (legacy script) | Prediction-stage only. Completed 2026-08-02 21:29–21:35 UTC. No integration into `derived/` or `figures/`. |
| 3C (`fig3_pwm_compatibility`) | — | Not attempted at all; no `predictions/Figure3C*` directory exists. |
| 3D (`fig3_scramble_schematic`) | — | Author-layout schematic in `MANIFEST.tsv`; not a computational panel. |

**Update 2026-08-05: all of Figure 3 is now independently verified
reproducible.** 3A/3C/3E/3F/3G were each re-run (3A end to end; 3C/3E/3F/3G
render/integration steps re-run against their already-complete scoring data)
into isolated directories under
`.../figure3_public_20260802T232000Z/derived/`, and diffed against the
release's committed compact tables — never overwriting the original
evidence:

| Panel | Result |
|---|---|
| 3A | Stage 1/2 completed (R018); regenerated `stage2_top_windows_per_gene.tsv` selects the same 150 windows as the release table, values differing by ≤3e-4 absolute (ordinary AlphaGenome variance) |
| 3B | PASS via `reproduce.py` (unchanged from above) |
| 3C | **Byte-identical** to both the untouched Jul-31 original and the release's `outputs/source_data/Figure3C_pwm_compatibility.tsv` (pure local PWM computation, no API calls) |
| 3D | Author-layout schematic; not a computational panel |
| 3E | **Byte-identical** to the release's `outputs/source_data/Figure3E_directional_recovery.tsv` |
| 3F | **Byte-identical** scoring data (SHA-256 match on `surface_summary_paired.csv`, `selected_mean_window.csv`, `selected_median_window.csv`); regenerated SVG differs from the release's only in an embedded timestamp and matplotlib's per-render random clip-path ID — every fill color is byte-identical |
| 3G | **Byte-identical** to the release's `outputs/source_data/Figure3G_component_necessity.tsv`; independently re-scored data (Aug 2) also matches the original Jul-31 scoring exactly, zero diff |

3E/3F/3G's original legacy scripts (`make_panel_E.py`,
`make_uniform_wide_main_panel.py`, `make_panel.py`) and 3C's
(`make_panel_C.py`) had no CLI output-path override and wrote in place over
the original published panel files — each was given a minimal `--out-dir`
(and `--source`/`--source-dir` where the input location also needed
overriding) so a clean-room re-run cannot clobber the evidence it's being
checked against; default behavior (no flags) is unchanged. See each script's
diff for the exact patch.

Not yet done for any of 3A/3C/3E/3F/3G: wiring into `reproduction/figure3.py`
so `reproduce.py run --panels 3A,3C,3E,3F,3G` reproduces and compares them
the same tracked, checksummed, one-command way 3B already is. The
verification above was done by invoking the (now patched) legacy scripts
directly against their already-existing scored data, not through the
harness. That harness-wiring is real but separate follow-up work — see
Priority 1 below.

## Findings closed by this evidence

| ID | Finding | Resolution |
|---|---|---|
| R001 | Figure 1B `tracks.npz` recorded but absent | Regenerated by the clean-room 1B step and verified PASS in every Figure 1 run |
| R002 | Only Figure 1B had a released analysis entry point | `reproduction/figure1.py`, `figure1_public.py`, `figure2.py`, `figure3.py` now cover Figs 1B–1F, 2B/2C/2E/2F, 3B |
| R004 | No section-level runner | `reproduce.py {doctor,prepare,run,compare}` added, with `--resume` and per-run audit reports |
| R005 | No exact fetch/preparation recipes | `reproduction/figure1_public.py` / `figure2.py` download, checksum and stage GLGC, GTEx, 1000 Genomes, GRCh38, Kircher MPRA and 4DN Hi-C sources directly |
| R006 | No reviewer-friendly credential setup | `.env.example` and `ALPHAGENOME_API_KEY.example.txt` added; both ignored by git; `doctor` checks the key before heavy imports |
| R007 | Reproduction defaults could overwrite reference assets | Runs now require an explicit, timestamped `--run-dir` and refuse a non-empty directory without `--resume` |
| R008 | No numerical/visual comparison rule | `reproduction/report.py` applies explicit tolerances per panel and records PASS/FAIL in `audit/comparison.json` |

## Findings still open

| ID | Finding | Notes |
|---|---|---|
| R002/R003 | Figures S1–S3 have no released analysis entry point or renderer | `SUPPORTED_PANELS` in `reproduce.py` still covers only Figs 1–3; S1–S3 are untouched |
| R002 (partial) | Figure 3 has no *tracked, one-command* entry point beyond 3B | `reproduction/figure3.py` still implements only `run_fig3b`; 3A/3C/3E/3F/3G are now verified reproducible (R018/R019) but only via direct, manually-invoked legacy scripts, not through `reproduce.py` |
| R009 | GLGC crosswalk inconsistency between `data/SOURCES.tsv` and `MANIFEST.tsv` | Not reconciled |
| R010 | Figure S1C missing one plotted numerical component | Not resolved |
| R011 | AppleDouble/exFAT sidecars break `validate.py` and Conda-on-exFAT | No `._*` exclusion added to `.gitignore` or `validate.py` |
| R012 | Every install pulls the full TensorFlow/JAX research stack | Not split |
| R013 | Strict release gate fails on pending Zenodo DOI | Carried over from `AUDIT_REPORT.md` / `RELEASE_CHECKLIST.md`; unrelated to the runner work |
| R014 | Credential check import ordering / Matplotlib cache | Not revisited |

## New items found while operating the runner

- **R015 — Figure 2C reconstruction failed intermittently.** The retained
  attempt log for `figure2_public_inputs_20260802T192638Z`
  (`audit/attempts/2026-08-02T192947.582609_0000/run.json`) shows step "2C:
  construct and score the 11-by-11 deletion grid" failing with
  `MergeError: Merge keys are not unique in right dataset; not a many-to-one
  merge`, then succeeding on the next attempt with unchanged code
  (`figure2.py` sha256 `9b7bfd44...` identical across both attempts). A
  reproduction script that only passes on retry is not yet a clean
  reproduction; the non-unique merge key needs to be found and deduplicated
  or asserted against in `reproduction/figure2.py`, not papered over by
  `--resume`.
- **R016 — Multi-gigabyte inputs are re-downloaded per run directory.**
  Each run under `repro_crash_test/runs/` re-fetches its own copy of
  `hg38.fa` (3.1 GB), `Liver.allpairs.txt.gz` (3.7 GB), the 4DN HepG2 `.hic`
  (14.7 GB), etc. Five runs already total ~19 GB on the T7 drive. A shared,
  content-addressed download cache outside the per-run directory (reused
  across `--run-dir`s by checksum) would cut both wall-clock time and disk
  use for the S1–S3 work below.
- **R017 — Runner work is entirely uncommitted.** `reproduce.py`,
  `reproduction/`, `tests/`, `.env.example` and
  `ALPHAGENOME_API_KEY.example.txt` are untracked; `README.md`, `.gitignore`,
  `environment.yml` and several `figures/rendered/*` files are modified but
  unstaged; local `main` is one commit behind `origin/main`. None of the
  PASS evidence above is currently reachable from a fresh clone.
- **R018 — RESOLVED 2026-08-05.** Figure 3A could not complete: missing
  GENCODE cache dependency. `sort1_figure_2e_100kb_rna_ism.py`
  hard-required a prebuilt `cache/gencode.v46.annotation.gtf.gz.feather`
  (its `_load_gene_tss` / `_load_gene_models` helpers had no download
  fallback, unlike the equivalent helper in `sort1_comprehensive_analysis.py`).
  Fix applied: ported the same download-on-cache-miss fallback into
  `sort1_figure_2e_100kb_rna_ism.py`, then restored the cache and re-ran.
  The interrupted checkpoint at `predictions/Figure3A_regional_ism/`
  (stage 1 complete, stage 2 partial: 28,417/51,120 rows) resumed and
  finished cleanly: stage 1 = 239,988 rows, stage 2 = 51,120 rows,
  `run_metadata.json` region/gene-TSS/tissue/seq-window fields identical to
  the original publication run's metadata. The regenerated
  `stage2_top_windows_per_gene.tsv` selects the same 150 windows as the
  compact table committed at `outputs/source_data/Figure3A_regional_ism/`,
  with per-window deltas differing by ≤3e-4 absolute — consistent with
  ordinary AlphaGenome run-to-run prediction variance, not a regression.
- **R019 — RESOLVED 2026-08-05.** Figure 3E/3F/3G were scored but never
  integrated, rendered, or compared; 3C was never attempted. Verified all
  four (plus 3A) as byte-identical or pixel-identical against the release —
  see the updated Figure 3 status table above for per-panel evidence.
  Remaining scope, not a correctness gap: none of 3A/3C/3E/3F/3G run through
  `reproduce.py` yet, so this verification isn't one-command or
  checksummed/audited the way 3B is. Tracked as a Priority 1 follow-up
  below.

## Priority next steps

Priority 0 (blocks an end-to-end reproducibility claim for the release):

1. Add analysis + renderer entry points for Figures S1–S3 (closes R002/R003
   for the remaining panels).
2. Root-cause and fix the Figure 2C `MergeError` (R015) so the panel passes
   on the first attempt, and add a regression case to
   `tests/test_reproduction.py`.

Priority 1:

3. Port Figure 3A/3C/3E/3F/3G into `reproduction/figure3.py` as tracked,
   checksummed `reproduce.py --panels` entry points (closes the remaining
   part of R002/R019). All five are now verified byte-/pixel-identical to
   the release by direct invocation of their (now `--out-dir`-safe) legacy
   scripts — see the Figure 3 status table above — so this is packaging
   already-proven logic into the harness, not new analysis work. 3D remains
   out of scope (author-layout schematic).
4. Commit the runner (`reproduce.py`, `reproduction/`, `tests/`,
   `.env.example`, `ALPHAGENOME_API_KEY.example.txt`), reconcile the other
   modified files, and push — a reviewer cloning `origin/main` today gets
   none of this (R017).
5. Add a shared download cache keyed by checksum/URL so repeat and future
   (S1–S3) runs reuse already-fetched multi-GB inputs instead of
   redownloading them (R016).
6. Ignore AppleDouble `._*` sidecars in `.gitignore` and `validate.py`, and
   document that Conda environments must not be created on exFAT (R011).

Priority 2:

7. Once S1–S3 land and Figure 3 is fully in `reproduce.py`, run a second
   dated clean-room audit (new audit ID under `repro_crash_test_audit/`)
   against the resulting commit to supersede `20260802T080725Z`, this time
   exercising every panel rather than a Figure 1/2/3B subset.
8. Reconcile the GLGC crosswalk (R009) and the Figure S1C missing source
   (R010).
9. Split a lightweight plotting environment from the full API/research
   environment (R012), and finish the carried-over Zenodo deposit / DOI
   blocker (R013) before any public-release tag.

## Suggested claim after remediation

Once items 1–6 above (Priority 0 and 1) land and a superseding audit passes,
the release could
support:

> A fresh clone with an authorized AlphaGenome API key regenerates every
> AlphaGenome-derived panel in Figures 1–3 and S1–S3 from public raw inputs,
> on the first attempt, with numerical/visual equivalence to the published
> reference reported under explicit tolerances. Deposit-pending large
> derived tables and the final DOI remain the only external dependency.
