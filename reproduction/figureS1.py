"""Clean-room AlphaGenome reconstruction for Supplementary Figure S1
(A-D): observed and predicted 3D contact architecture around rs12740374.

Ported from the working archive's report/panel_contact/
run_observed_hic_validation.py, run_contact_architecture_analysis.py, and
results/figure1/supplementary/contact_maps_rs12740374/
run_figure1_supplementary_contact_maps_rs12740374.py -- none of which are
part of this repository.

S1A and S1B need no new downloads or AlphaGenome calls: they are
re-derived from the raw Hi-C/FOLD_0 contact matrices Figure 1D/1E already
fetch and cache (predictions/Figure1D_observed_hic/observed.npz,
predictions/Figure1E_fold0_contact/ref_alt.npz), reusing the exact
same-distance-null statistic as the archive (verified against the
archive's own helper functions rather than assumed equivalent to
Figure 1F's superficially-similar percentile table, which uses a
different null-window formula).

S1C and S1D are real, fresh AlphaGenome calls: S1C is one official
ContactMapScorer call, one HepG2 ALL_FOLDS WT/ALT virtual-4C delta, and
100 geometry-matched null SNVs (deterministically sampled, seed
12740374, same RNG call sequence as the archive); S1D is one reference-
only ALL_FOLDS CONTACT_MAPS call across every available ontology.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .common import Audit, api_key
from .figure1 import GENES, HEPG2, VARIANT_ALT, VARIANT_CHR, VARIANT_POS, VARIANT_REF, VARIANT_RSID, _ag, _save_svg
from .figure1_public import run_observed_contact
from .figure4 import fetch_ucsc_hg38

TSS = {"CELSR2": 109_249_538, "PSRC1": 109_283_186, "SORT1": 109_397_918}
GENE_STRAND = {"CELSR2": "+", "PSRC1": "-", "SORT1": "-"}
GENE_BOUNDS = {
    "CELSR2": (109_249_538, 109_275_751),
    "PSRC1": (109_279_555, 109_283_186),
    "SORT1": (109_309_567, 109_397_918),
}
NULL_WINDOW_BP = 2 ** 20
NULL_ANCHOR_NEIGHBORHOOD_BP = 10_000
NULL_N_SAMPLES = 100
NULL_SEED = 12_740_374
NULL_PROMOTER_GENE = "SORT1"


# --- shared statistics helpers, ported verbatim from the archive's
# run_contact_architecture_analysis.py ------------------------------------

def _bin_for_1based(position: int, interval_start0: int, resolution: int) -> int:
    return int((int(position) - 1 - int(interval_start0)) // int(resolution))


def _virtual4c_row(matrix: np.ndarray, anchor_bin: int) -> np.ndarray:
    return 0.5 * (matrix[anchor_bin, :] + matrix[:, anchor_bin])


def _window_mean(vector: np.ndarray, target_bin: int, flank_bins: int = 1) -> float:
    lo = max(0, int(target_bin) - int(flank_bins))
    hi = min(vector.size, int(target_bin) + int(flank_bins) + 1)
    return float(np.nanmean(vector[lo:hi]))


def _same_distance_null(matrix: np.ndarray, separation_bins: int, *, flank_bins: int = 1) -> np.ndarray:
    d = abs(int(separation_bins))
    n = int(matrix.shape[0])
    values: list[float] = []
    for anchor in range(0, n - d):
        target = anchor + d
        if target - flank_bins < 0 or target + flank_bins >= n:
            continue
        values.append(_window_mean(_virtual4c_row(matrix, anchor), target, flank_bins))
    return np.asarray(values, dtype=float)


def _percentile_midrank(values: np.ndarray, observed: float) -> float:
    less = np.count_nonzero(values < observed)
    equal = np.count_nonzero(values == observed)
    return 100.0 * float(less + 0.5 * equal) / float(values.size)


def _promoter_rows(matrix: np.ndarray, *, interval_start0: int, resolution: int, anchor_bin: int, source_label: str) -> pd.DataFrame:
    v4c = _virtual4c_row(matrix, anchor_bin)
    rows: list[dict[str, Any]] = []
    for gene in GENES:
        tss = TSS[gene]
        target_bin = _bin_for_1based(tss, interval_start0, resolution)
        observed = _window_mean(v4c, target_bin, 1)
        null = _same_distance_null(matrix, target_bin - anchor_bin, flank_bins=1)
        null_sd = float(np.std(null, ddof=1))
        rows.append({
            "source": source_label, "gene": gene, "strand": GENE_STRAND[gene],
            "gene_start_hg38": GENE_BOUNDS[gene][0], "gene_end_hg38": GENE_BOUNDS[gene][1],
            "tss_hg38": tss, "tss_bin": target_bin,
            "tss_bin_start_0based": interval_start0 + target_bin * resolution,
            "tss_bin_end_0based": interval_start0 + (target_bin + 1) * resolution,
            "variant_to_tss_bp": int(tss - VARIANT_POS), "separation_bins": int(target_bin - anchor_bin),
            "log_observed_expected_pm1bin": observed, "implied_observed_expected": float(np.exp(observed)),
            "same_distance_null_n": int(null.size), "same_distance_null_mean": float(np.mean(null)),
            "same_distance_null_sd": null_sd, "same_distance_z": float((observed - np.mean(null)) / null_sd),
            "same_distance_percentile": _percentile_midrank(null, observed),
            "same_distance_high_empirical_p": float((np.count_nonzero(null >= observed) + 1) / (null.size + 1)),
        })
    return pd.DataFrame(rows)


# --- S1A: observed-only virtual-4C track (zero new fetches) ---------------

def run_figs1a(run_dir: Path, audit: Audit) -> None:
    from .figure1 import run_fig1e

    cache = run_dir / "predictions/Figure1D_observed_hic/observed.npz"
    fold0_cache = run_dir / "predictions/Figure1E_fold0_contact/ref_alt.npz"
    if not cache.exists():
        # run_observed_contact (Figure 1D) also builds Figure 1F internally,
        # which reads Figure 1E's cache -- that must exist first.
        if not fold0_cache.exists():
            with audit.step("S1A prerequisite: Figure 1E FOLD_0 contact map"):
                run_fig1e(run_dir, audit)
        with audit.step("S1A prerequisite: Figure 1D observed Hi-C"):
            run_observed_contact(run_dir, audit, None)
    with audit.step("S1A: observed SNP-anchored virtual-4C track"):
        z = np.load(cache)
        matrix, start, resolution = np.asarray(z["matrix"], dtype=float), int(z["start"]), int(z["resolution"])
        anchor = _bin_for_1based(VARIANT_POS, start, resolution)
        v4c = _virtual4c_row(matrix, anchor)
        # anchor_flank=1: average matrix rows/cols anchor-1..anchor+1 (not just the single anchor bin).
        v4c_flanked = 0.5 * (
            np.nanmean(matrix[max(0, anchor - 1): anchor + 2, :], axis=0)
            + np.nanmean(matrix[:, max(0, anchor - 1): anchor + 2], axis=1)
        )
        smoothed = np.convolve(v4c_flanked, np.ones(3) / 3.0, mode="same")
        centers = start + (np.arange(matrix.shape[0]) + 0.5) * resolution
        offset_kb = (centers - (VARIANT_POS - 1)) / 1000.0
        keep = (offset_kb >= -55) & (offset_kb <= 165)
        out = run_dir / "derived/FigureS1A_observed_virtual4c.tsv"
        out.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({
            "position_from_rs12740374_kb": offset_kb[keep],
            "observed_contact_oe_3x3_smoothed": smoothed[keep],
        }).to_csv(out, sep="\t", index=False)
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(4, 2.2))
        ax.plot(offset_kb[keep], smoothed[keep], color="#3F6FA8", lw=1.0)
        ax.axvline(0, color="#222", ls=":", lw=0.8)
        ax.set_xlabel("Position from rs12740374 (kb)"); ax.set_ylabel("Observed O/E (3x3 smoothed)")
        _save_svg(fig, run_dir / "figures/FigureS1A.svg")


# --- S1B: observed + FOLD_0 TSS-bin contact, same-distance percentile -----

def _contact_window(matrix: np.ndarray, anchor: int, target: int, flank: int) -> float:
    return float(np.mean(matrix[anchor - flank: anchor + flank + 1, target - flank: target + flank + 1]))


def _same_distance_percentile_raw(matrix: np.ndarray, anchor: int, target: int, flank: int) -> tuple[float, float, int]:
    """S1B's actual statistic (build_figure1_supplementary_panels.py::
    contact_resolution_sensitivity) -- a raw single-bin (flank=0) value and
    a same-distance null over ungapped left/right bin pairs starting at 0,
    not the virtual-4C-row statistic S1C/S1D use. Verified byte-identical
    against the archive's own frozen S1A output before trusting this."""
    observed_value = _contact_window(matrix, anchor, target, flank)
    distance = abs(target - anchor)
    null = np.asarray([_contact_window(matrix, left, left + distance, flank) for left in range(flank, matrix.shape[0] - distance - flank)])
    return observed_value, _percentile_midrank(null, observed_value), int(null.size)


def run_figs1b(run_dir: Path, audit: Audit) -> None:
    from .figure1 import run_fig1e

    obs_cache = run_dir / "predictions/Figure1D_observed_hic/observed.npz"
    fold0_cache = run_dir / "predictions/Figure1E_fold0_contact/ref_alt.npz"
    if not obs_cache.exists():
        with audit.step("S1B prerequisite: Figure 1D observed Hi-C"):
            run_observed_contact(run_dir, audit, None)
    if not fold0_cache.exists():
        with audit.step("S1B prerequisite: Figure 1E FOLD_0 contact map"):
            run_fig1e(run_dir, audit)
    with audit.step("S1B: TSS-bin contact, observed vs. AlphaGenome FOLD_0"):
        zo = np.load(obs_cache)
        obs_matrix, obs_start, obs_res = np.asarray(zo["matrix"], dtype=float), int(zo["start"]), int(zo["resolution"])
        zf = np.load(fold0_cache)
        fold0_matrix, fold0_start, fold0_res = np.asarray(zf["ref"], dtype=float), int(zf["start"]), int(zf["resolution"])

        rows = []
        for label, matrix, start, resolution in (
            ("Experimental HepG2 Hi-C", obs_matrix, obs_start, obs_res),
            ("AlphaGenome FOLD_0", fold0_matrix, fold0_start, fold0_res),
        ):
            anchor = _bin_for_1based(VARIANT_POS, start, resolution)
            for gene in GENES:
                target = _bin_for_1based(TSS[gene], start, resolution)
                value, percentile, null_n = _same_distance_percentile_raw(matrix, anchor, target, flank=0)
                rows.append({
                    "source": label, "summary": "TSS-containing bin", "gene": gene, "resolution_bp": resolution,
                    "contact_value": value, "same_distance_percentile": percentile, "same_distance_null_n": null_n,
                })
        combined = pd.DataFrame(rows)
        out = run_dir / "derived/FigureS1B_tss_bin_contact.tsv"
        out.parent.mkdir(parents=True, exist_ok=True)
        combined.to_csv(out, sep="\t", index=False)


# --- S1C: rs12740374 allele-specific contact delta + null -----------------

def _sample_null_snvs(fasta, *, chrom: str, center_pos1: int, half_width_bp: int, exclude_pos1: int, n: int, seed: int) -> list[tuple[int, str, str]]:
    start0, end0 = center_pos1 - 1 - half_width_bp, center_pos1 - 1 + half_width_bp
    seq = fasta.fetch(chrom, start0, end0).upper()
    rng = np.random.default_rng(int(seed))
    valid_idx = np.array([i for i, b in enumerate(seq) if b in "ACGT" and (start0 + i + 1) != exclude_pos1], dtype=int)
    picked = rng.choice(valid_idx, size=min(n, valid_idx.size), replace=False)
    bases = np.array(list("ACGT"), dtype="<U1")
    out: list[tuple[int, str, str]] = []
    for i in picked:
        ref = str(seq[int(i)])
        alt = str(rng.choice(bases[bases != ref]))
        out.append((int(start0 + int(i) + 1), ref, alt))
    return out


def run_figs1c(run_dir: Path, audit: Audit) -> None:
    genome, ontology, dna_client, model_modules = _ag()
    dna_model, variant_scorers = model_modules
    out = run_dir / "derived/FigureS1C_contact_allele_delta"
    out.mkdir(parents=True, exist_ok=True)
    interval = genome.Interval(VARIANT_CHR, VARIANT_POS, VARIANT_POS).resize(NULL_WINDOW_BP)

    # (a) official multi-track ContactMapScorer for rs12740374.
    scores_path = out / "panelC_contact_map_scores.tsv"
    if not scores_path.exists():
        with audit.step("S1C: official ContactMapScorer for rs12740374"):
            client = dna_client.create(api_key(), model_version=dna_model.ModelVersion.ALL_FOLDS, timeout=300)
            variant = genome.Variant(VARIANT_CHR, VARIANT_POS, VARIANT_REF, VARIANT_ALT, VARIANT_RSID)
            scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["CONTACT_MAPS"]
            scores = client.score_variant(interval, variant, [scorer])
            df = variant_scorers.tidy_scores(scores)
            df.to_csv(scores_path, sep="\t", index=False)
            audit.add_api_calls("S1C", 1)
            audit.add_api_requests("S1C", 1)

    # (b) HepG2 WT/ALT virtual-4C promoter deltas.
    delta_path = out / "panelB_promoter_allele_delta.tsv"
    hepg2_cache = run_dir / "predictions/FigureS1C_hepg2_wt_alt/ref_alt.npz"
    if not hepg2_cache.exists():
        with audit.step("S1C: HepG2 ALL_FOLDS WT/ALT contact map for rs12740374"):
            client = dna_client.create(api_key(), model_version=dna_model.ModelVersion.ALL_FOLDS, timeout=300)
            variant = genome.Variant(VARIANT_CHR, VARIANT_POS, VARIANT_REF, VARIANT_ALT, VARIANT_RSID)
            output = client.predict_variant(
                interval=interval, variant=variant, requested_outputs=[dna_client.OutputType.CONTACT_MAPS],
                ontology_terms=[ontology.from_curie(HEPG2)],
            )
            ref = np.nanmean(np.asarray(output.reference.contact_maps.values), axis=-1).astype(np.float32)
            alt = np.nanmean(np.asarray(output.alternate.contact_maps.values), axis=-1).astype(np.float32)
            hepg2_cache.parent.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(hepg2_cache, ref=ref, alt=alt, start=int(interval.start), resolution=int(output.reference.contact_maps.resolution))
            audit.add_api_calls("S1C", 1)
            audit.add_api_requests("S1C", 1)
    if not delta_path.exists():
        z = np.load(hepg2_cache)
        ref, alt = np.asarray(z["ref"], dtype=float), np.asarray(z["alt"], dtype=float)
        start, resolution = int(z["start"]), int(z["resolution"])
        delta_matrix = alt - ref
        anchor = _bin_for_1based(VARIANT_POS, start, resolution)
        v4c_delta = _virtual4c_row(delta_matrix, anchor)
        rows = []
        for gene in GENES:
            target_bin = _bin_for_1based(TSS[gene], start, resolution)
            delta = _window_mean(v4c_delta, target_bin, 1)
            rows.append({"gene": gene, "delta_log_observed_expected_pm1bin": delta, "percent_change_implied_observed_expected": 100.0 * (np.exp(delta) - 1.0)})
        pd.DataFrame(rows).to_csv(delta_path, sep="\t", index=False)

    # (c) 100 geometry-matched null SNVs around rs12740374.
    null_path = out / "panelC_local_snv_null.tsv"
    if not null_path.exists():
        z = np.load(hepg2_cache)
        start, resolution = int(z["start"]), int(z["resolution"])
        promoter_bin = _bin_for_1based(TSS[NULL_PROMOTER_GENE], start, resolution)
        z_delta = np.asarray(z["alt"], dtype=float) - np.asarray(z["ref"], dtype=float)
        anchor0 = _bin_for_1based(VARIANT_POS, start, resolution)
        observed_metric = _window_mean(_virtual4c_row(z_delta, anchor0), promoter_bin, 1)

        import pysam
        ucsc_fasta_path = fetch_ucsc_hg38(run_dir, audit)
        fasta = pysam.FastaFile(str(ucsc_fasta_path))
        samples = _sample_null_snvs(
            fasta, chrom=VARIANT_CHR, center_pos1=VARIANT_POS, half_width_bp=NULL_ANCHOR_NEIGHBORHOOD_BP,
            exclude_pos1=VARIANT_POS, n=NULL_N_SAMPLES, seed=NULL_SEED,
        )
        client = dna_client.create(api_key(), model_version=dna_model.ModelVersion.ALL_FOLDS, timeout=300)
        rows = []
        with audit.step(f"S1C: score {len(samples)} null SNVs for the geometry-matched contact-delta distribution"):
            for pos1, ref_base, alt_base in samples:
                var = genome.Variant(VARIANT_CHR, pos1, ref_base, alt_base, f"null_{pos1}_{ref_base}>{alt_base}")
                output = client.predict_variant(
                    interval=interval, variant=var, requested_outputs=[dna_client.OutputType.CONTACT_MAPS],
                    ontology_terms=[ontology.from_curie(HEPG2)],
                )
                ref_arr = np.nanmean(np.asarray(output.reference.contact_maps.values), axis=-1)
                alt_arr = np.nanmean(np.asarray(output.alternate.contact_maps.values), axis=-1)
                delta_arr = alt_arr - ref_arr
                anchor_bin = _bin_for_1based(pos1, start, resolution)
                delta_v4c = _virtual4c_row(delta_arr, anchor_bin)
                metric = _window_mean(delta_v4c, promoter_bin, 1)
                rows.append({
                    "variant_pos_hg38": pos1, "ref": ref_base, "alt": alt_base, "anchor_bin_global": anchor_bin,
                    "promoter_gene": NULL_PROMOTER_GENE, "promoter_bin_global": promoter_bin,
                    "delta_contact_virtual4c_pm1bin": metric, "ontology_curie": HEPG2, "window_bp": NULL_WINDOW_BP,
                    "anchor_neighborhood_bp": NULL_ANCHOR_NEIGHBORHOOD_BP, "anchor_distance_from_rs_bp": pos1 - VARIANT_POS,
                    "seed": NULL_SEED, "observed_rs12740374_delta_contact_virtual4c_pm1bin": observed_metric,
                })
            audit.add_api_calls("S1C", len(samples))
            audit.add_api_requests("S1C", len(samples))
        pd.DataFrame(rows).to_csv(null_path, sep="\t", index=False)


# --- S1D: reference-only contact map across every ontology ----------------

def run_figs1d(run_dir: Path, audit: Audit) -> None:
    genome, ontology, dna_client, model_modules = _ag()
    dna_model, _ = model_modules
    cache = run_dir / "predictions/FigureS1D_all_contexts/values.npz"
    if not cache.exists():
        with audit.step("S1D: reference-only ALL_FOLDS contact map across every ontology"):
            client = dna_client.create(api_key(), model_version=dna_model.ModelVersion.ALL_FOLDS, timeout=300)
            interval = genome.Interval(VARIANT_CHR, VARIANT_POS, VARIANT_POS).resize(NULL_WINDOW_BP)
            output = client.predict_interval(interval=interval, requested_outputs={dna_client.OutputType.CONTACT_MAPS}, ontology_terms=None)
            tracks = output.contact_maps
            values = np.asarray(tracks.values, dtype=np.float32)
            cache.parent.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(cache, values=values, start=int(tracks.interval.start), resolution=int(tracks.resolution))
            tracks.metadata.reset_index(drop=True).to_csv(cache.parent / "track_metadata.tsv", sep="\t", index=False)
            audit.add_api_calls("S1D", 1)
            audit.add_api_requests("S1D", 1)
    z = np.load(cache)
    values, start, resolution = np.asarray(z["values"], dtype=float), int(z["start"]), int(z["resolution"])
    metadata = pd.read_csv(cache.parent / "track_metadata.tsv", sep="\t")
    anchor = _bin_for_1based(VARIANT_POS, start, resolution)
    context_rows = []
    for (ontology_curie, biosample_name), idx in metadata.groupby(["ontology_curie", "biosample_name"], sort=False).groups.items():
        indices = list(idx)
        mean_map = np.nanmean(values[:, :, indices], axis=2)
        promoter = _promoter_rows(mean_map, interval_start0=start, resolution=resolution, anchor_bin=anchor, source_label=f"{biosample_name} mean")
        promoter.insert(1, "ontology_curie", str(ontology_curie))
        promoter.insert(2, "biosample_name", str(biosample_name))
        promoter.insert(3, "n_tracks", len(indices))
        context_rows.append(promoter)
    context_df = pd.concat(context_rows, ignore_index=True)
    out = run_dir / "derived/FigureS1D_contact_contexts.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    context_df.to_csv(out, sep="\t", index=False)
