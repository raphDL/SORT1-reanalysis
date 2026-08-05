"""Clean-room AlphaGenome reconstruction for Figure 4E/4F: chromosome-1
Hi-C-guided distal-enhancer transfer benchmark (real observed HepG2 Hi-C
contact strength vs. AlphaGenome-predicted T/G transfer effect).

Ported from the working archive's report/panel_distal_hic_transfer/
(build_chr1_promoter_hic_catalog.py, run_chr1_distal_315_transfer.py,
analyze_chr1_distal_transfer.py, plot_distal_contact_results.py) -- none of
which are part of this repository. See REPRODUCIBILITY_NEXT_STEPS.md.

Two independent stages:

1. Promoter catalogue + high/low contact-site selection (this file's
   `select_sites_*` machinery). This touches GENCODE and an observed 4DN
   HepG2 Hi-C map -- no AlphaGenome calls at all -- and is fully
   deterministic given those two public, versioned sources. Re-derived at
   run time every time, exactly like Figure 4C's HPA cohort selection (not
   frozen -- the working archive's own build script is directly portable).

2. AlphaGenome transfer scoring (`_score_distal_states` /
   `run_fig4ef`): one native sequence plus minor-T/major-G 315bp
   replacements at each of the 3 selected sites per promoter (7 states),
   scored for RNA_SEQ + ATAC + CHIP_HISTONE:H3K27ac in HepG2. Every
   prediction is freshly computed (checkpointed per sequence hash so an
   interrupted run can resume, but never seeded from previously-scored
   data).

The strict per-modality track selector (reproduction/data/
figure4ef_track_selection.json) is a frozen *design* input (which exact
HepG2 tracks count as "the" RNA/ATAC/H3K27ac tracks), not an AlphaGenome
prediction -- same category as the JASPAR PFM file or Figure 4B's frozen
recipient list.
"""

from __future__ import annotations

import hashlib
import json
import math
from itertools import combinations
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd

from .common import Audit, api_key, download
from .figure1 import _save_svg
from .figure2 import HEPG2
from .figure4 import DATA_DIR, GENCODE_URL, _ensure_gencode, fetch_ucsc_hg38, make_asymmetric_315_donors

CHROM = "chr1"
HIC_CHROM = "1"
HIC_URL = (
    "https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/"
    "wfoutput/25104375-a588-46e6-a382-663cee6c332f/4DNFICSTCJQZ.hic"
)
TRACK_SPEC_PATH = DATA_DIR / "figure4ef_track_selection.json"

PROMOTER_HALF_WIDTH = 5_000
MIN_DISTANCE = 50_000
MAX_DISTANCE = 1_000_000
RESOLUTION = 5_000
DONOR_LENGTH = 315
DONOR_VARIANT_INDEX = 179
MODEL_LENGTH = 2 ** 20
LOW_CONTROLS = 2
LOW_DISTANCE_CALIPER = 75_000
MIN_SITE_SEPARATION = 15_000
LOW_GC_TOLERANCE = 0.08
HIGH_MIN_OE = 1.0
HIGH_MIN_PERCENTILE = 0.90
DISTANCE_BAND_WIDTH = 50_000
LOCAL_CONTEXT_LENGTH = 5_000
TSS_HALF_WIDTH = 2_000
LOCAL_HALF_WIDTH = 1_000
SITE_ROLES = ("high_contact", "low_contact_1", "low_contact_2")
EDIT_CONSTRUCTS = ("minor_T", "major_G")
PRIMARY_LOW_ROLES = ("low_contact_1", "low_contact_2")
PRIMARY_MAX_DISTANCE = 500_000
PRIMARY_MAX_N_FRACTION = 0.01
BOOTSTRAP_SEED = 12_740_374
BOOTSTRAP_N = 2_000


# --- stage 1: chr1 promoter catalogue + observed Hi-C site selection ------

def _load_promoters(gencode_path: Path, *, chromosome_length: int) -> pd.DataFrame:
    cols = ["Chromosome", "Feature", "Start", "End", "Strand", "gene_name", "gene_id", "gene_type"]
    frame = pd.read_feather(gencode_path, columns=cols)
    frame = frame[frame["Feature"].astype(str).eq("gene") & frame["Chromosome"].astype(str).eq(CHROM)].copy()
    frame["start0"] = pd.to_numeric(frame["Start"], errors="raise").astype(int)
    frame["end0"] = pd.to_numeric(frame["End"], errors="raise").astype(int)
    frame["tss0"] = np.where(frame["Strand"].astype(str).eq("+"), frame["start0"], frame["end0"] - 1).astype(int)
    frame["tss1"] = frame["tss0"] + 1
    frame["promoter_start0"] = (frame["tss0"] - PROMOTER_HALF_WIDTH).clip(lower=0)
    frame["promoter_end0"] = (frame["tss0"] + PROMOTER_HALF_WIDTH).clip(upper=int(chromosome_length))
    frame["gene_name"] = frame["gene_name"].fillna("").astype(str)
    frame["gene_id"] = frame["gene_id"].fillna("").astype(str)
    frame["gene_id_base"] = frame["gene_id"].str.replace(r"\.[0-9]+$", "", regex=True)
    frame["benchmark_gene_type"] = frame["gene_type"].astype(str).eq("protein_coding")
    frame["promoter_id"] = frame["gene_id_base"] + "__" + frame["gene_name"] + "__tss" + frame["tss1"].astype(str)
    frame = frame.sort_values(["tss0", "gene_id", "gene_name"], kind="stable").reset_index(drop=True)
    frame.insert(0, "promoter_index", np.arange(len(frame), dtype=int))
    return frame.rename(columns={"Chromosome": "chrom", "Strand": "strand"})[
        ["promoter_index", "promoter_id", "chrom", "gene_id", "gene_id_base", "gene_name", "gene_type",
         "strand", "start0", "end0", "tss0", "tss1", "promoter_start0", "promoter_end0", "benchmark_gene_type"]
    ]


def _load_transcript_promoters(gencode_path: Path, *, chromosome_length: int) -> pd.DataFrame:
    cols = ["Chromosome", "Feature", "Start", "End", "Strand"]
    frame = pd.read_feather(gencode_path, columns=cols)
    frame = frame[frame["Feature"].astype(str).eq("transcript") & frame["Chromosome"].astype(str).eq(CHROM)].copy()
    frame["start0"] = pd.to_numeric(frame["Start"], errors="raise").astype(int)
    frame["end0"] = pd.to_numeric(frame["End"], errors="raise").astype(int)
    # Assign the bare ndarray (positional) before deriving Series from it --
    # `frame` keeps its original, non-contiguous post-filter index, so
    # wrapping a fresh-0-based-index pd.Series() straight into a column
    # assignment here would silently index-align instead of assigning
    # positionally and corrupt most rows.
    frame["tss0"] = np.where(frame["Strand"].astype(str).eq("+"), frame["start0"], frame["end0"] - 1).astype(int)
    frame["promoter_start0"] = (frame["tss0"] - PROMOTER_HALF_WIDTH).clip(lower=0)
    frame["promoter_end0"] = (frame["tss0"] + PROMOTER_HALF_WIDTH).clip(upper=int(chromosome_length))
    return frame[["promoter_start0", "promoter_end0"]]


def _load_merged_exons(gencode_path: Path) -> tuple[np.ndarray, np.ndarray]:
    frame = pd.read_feather(gencode_path, columns=["Chromosome", "Feature", "Start", "End"])
    frame = frame[frame["Chromosome"].astype(str).eq(CHROM) & frame["Feature"].astype(str).eq("exon")].copy()
    return _merge_intervals(
        pd.to_numeric(frame["Start"], errors="raise").to_numpy(int),
        pd.to_numeric(frame["End"], errors="raise").to_numpy(int),
    )


def _merge_intervals(starts: np.ndarray, ends: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    order = np.argsort(starts, kind="stable")
    starts = np.asarray(starts, dtype=int)[order]
    ends = np.asarray(ends, dtype=int)[order]
    merged_starts: list[int] = []
    merged_ends: list[int] = []
    for start, end in zip(starts, ends, strict=True):
        if not merged_starts or int(start) > merged_ends[-1]:
            merged_starts.append(int(start))
            merged_ends.append(int(end))
        else:
            merged_ends[-1] = max(merged_ends[-1], int(end))
    return np.asarray(merged_starts, dtype=int), np.asarray(merged_ends, dtype=int)


def _overlaps_merged(starts: np.ndarray, ends: np.ndarray, merged_starts: np.ndarray, merged_ends: np.ndarray) -> np.ndarray:
    indices = np.searchsorted(merged_starts, ends - 1, side="right") - 1
    output = np.zeros(len(starts), dtype=bool)
    valid = indices >= 0
    output[valid] = merged_ends[indices[valid]] > starts[valid]
    return output


def _aligned_bins(start0: int, end0: int) -> np.ndarray:
    first = (int(start0) // RESOLUTION) * RESOLUTION
    last = ((int(end0) - 1) // RESOLUTION) * RESOLUTION
    return np.arange(first, last + RESOLUTION, RESOLUTION, dtype=int)


def _extract_promoter_contacts(
    zoom: Any, promoter: pd.Series, *, chromosome_length: int,
    merged_promoter_starts: np.ndarray, merged_promoter_ends: np.ndarray,
    merged_exon_starts: np.ndarray, merged_exon_ends: np.ndarray,
) -> pd.DataFrame:
    tss0 = int(promoter["tss0"])
    anchor_bins = np.asarray([(tss0 // RESOLUTION) * RESOLUTION], dtype=int)
    query_start0 = max(0, tss0 - MAX_DISTANCE - RESOLUTION)
    query_end0 = min(int(chromosome_length), tss0 + MAX_DISTANCE + RESOLUTION + 1)
    target_bins = _aligned_bins(query_start0, query_end0)
    target_index = {int(v): i for i, v in enumerate(target_bins)}
    anchor_set = {int(v) for v in anchor_bins}
    contact_sum = np.zeros(len(target_bins), dtype=float)
    observed_anchor_bins = np.zeros(len(target_bins), dtype=int)

    records = zoom.getRecords(
        int(anchor_bins[0]), int(anchor_bins[-1] + RESOLUTION - 1),
        int(target_bins[0]), int(target_bins[-1] + RESOLUTION - 1),
    )
    for record in records:
        bin_x, bin_y, value = int(record.binX), int(record.binY), float(record.counts)
        if not np.isfinite(value):
            continue
        if bin_x in anchor_set and bin_y in target_index and bin_y not in anchor_set:
            index = target_index[bin_y]
        elif bin_y in anchor_set and bin_x in target_index and bin_x not in anchor_set:
            index = target_index[bin_x]
        else:
            continue
        contact_sum[index] += value
        observed_anchor_bins[index] += 1

    contact = contact_sum / float(len(anchor_bins))
    smoothed = np.convolve(contact, np.ones(3, dtype=float) / 3.0, mode="same")
    target_centers0 = target_bins + RESOLUTION // 2
    signed_distance = target_centers0 - tss0
    absolute_distance = np.abs(signed_distance)
    keep = (
        (absolute_distance >= MIN_DISTANCE) & (absolute_distance <= MAX_DISTANCE)
        & (target_bins >= 0) & (target_bins + RESOLUTION <= int(chromosome_length))
    )
    target_bins, target_centers0 = target_bins[keep], target_centers0[keep]
    signed_distance, absolute_distance = signed_distance[keep], absolute_distance[keep]
    contact, smoothed, observed_anchor_bins = contact[keep], smoothed[keep], observed_anchor_bins[keep]
    target_ends = target_bins + RESOLUTION
    promoter_overlap = _overlaps_merged(target_bins, target_ends, merged_promoter_starts, merged_promoter_ends)
    replacement_start0 = target_centers0 - DONOR_LENGTH // 2
    replacement_end0 = replacement_start0 + DONOR_LENGTH
    exon_overlap = _overlaps_merged(replacement_start0, replacement_end0, merged_exon_starts, merged_exon_ends)
    target_gene_overlap = (replacement_start0 < int(promoter["end0"])) & (replacement_end0 > int(promoter["start0"]))
    frame = pd.DataFrame({
        "promoter_index": int(promoter["promoter_index"]), "promoter_id": str(promoter["promoter_id"]),
        "gene_id": str(promoter["gene_id"]), "gene_id_base": str(promoter["gene_id_base"]),
        "gene_name": str(promoter["gene_name"]), "gene_type": str(promoter["gene_type"]),
        "strand": str(promoter["strand"]), "tss0": tss0, "tss1": int(promoter["tss1"]),
        "target_bin_start0": target_bins, "target_bin_end0": target_ends, "target_bin_center0": target_centers0,
        "replacement_start0": replacement_start0, "replacement_end0": replacement_end0,
        "signed_distance_bp": signed_distance, "absolute_distance_bp": absolute_distance,
        "contact_oe": contact, "contact_oe_smooth3": smoothed, "n_nonzero_promoter_bins": observed_anchor_bins,
        "overlaps_any_10kb_promoter": promoter_overlap, "overlaps_any_exon": exon_overlap,
        "overlaps_target_gene_body": target_gene_overlap,
    })
    frame["eligible_nonpromoter_target"] = (
        ~frame["overlaps_any_10kb_promoter"] & ~frame["overlaps_any_exon"] & ~frame["overlaps_target_gene_body"]
        & frame["replacement_start0"].ge(0) & frame["replacement_end0"].le(int(chromosome_length))
        & np.isfinite(frame["contact_oe_smooth3"])
    )
    return frame


def _site_sequence_context(fasta, *, center0: int, chromosome_length: int) -> tuple[float, float, float, float]:
    start0 = int(center0) - LOCAL_CONTEXT_LENGTH // 2
    end0 = start0 + LOCAL_CONTEXT_LENGTH
    if start0 < 0 or end0 > int(chromosome_length):
        return math.nan, math.nan, math.nan, math.nan
    sequence = fasta.fetch(CHROM, start0, end0).upper()
    if len(sequence) != LOCAL_CONTEXT_LENGTH:
        return math.nan, math.nan, math.nan, math.nan
    slot_start = int(center0) - DONOR_LENGTH // 2 - start0
    slot_end = slot_start + DONOR_LENGTH
    slot, flank = sequence[slot_start:slot_end], sequence[:slot_start] + sequence[slot_end:]
    if len(slot) != DONOR_LENGTH or len(flank) != LOCAL_CONTEXT_LENGTH - DONOR_LENGTH:
        return math.nan, math.nan, math.nan, math.nan
    slot_canonical = sum(b in "ACGT" for b in slot)
    flank_canonical = sum(b in "ACGT" for b in flank)
    slot_gc = float(sum(b in "GC" for b in slot) / len(slot)) if slot_canonical == len(slot) and flank_canonical else math.nan
    flank_ambiguous = float(1.0 - flank_canonical / len(flank))
    canonical_flank = "".join(b for b in flank if b in "ACGT")
    if not canonical_flank:
        return slot_gc, math.nan, math.nan, flank_ambiguous
    flank_gc = float(sum(b in "GC" for b in canonical_flank) / len(canonical_flank))
    flank_cpg = float(canonical_flank.count("CG") / max(1, len(canonical_flank) - 1))
    return slot_gc, flank_gc, flank_cpg, flank_ambiguous


def _annotate_sequence_context(contacts: pd.DataFrame, *, fasta, chromosome_length: int) -> pd.DataFrame:
    out = contacts.copy()
    for col in ("native_replacement_gc", "local_flank_gc", "local_flank_cpg_density", "local_flank_ambiguous_fraction"):
        out[col] = math.nan
    base_mask = out["eligible_nonpromoter_target"].astype(bool)
    base_indices = out.index[base_mask].tolist()
    unique_centers = sorted(out.loc[base_mask, "target_bin_center0"].astype(int).unique().tolist())
    context_by_center = {
        int(c): _site_sequence_context(fasta, center0=int(c), chromosome_length=int(chromosome_length))
        for c in unique_centers
    }
    contexts = [context_by_center[int(out.at[i, "target_bin_center0"])] for i in base_indices]
    if base_indices:
        out.loc[base_indices, "native_replacement_gc"] = [v[0] for v in contexts]
        out.loc[base_indices, "local_flank_gc"] = [v[1] for v in contexts]
        out.loc[base_indices, "local_flank_cpg_density"] = [v[2] for v in contexts]
        out.loc[base_indices, "local_flank_ambiguous_fraction"] = [v[3] for v in contexts]
    out["eligible_sequence_target"] = (
        out["eligible_nonpromoter_target"].astype(bool)
        & np.isfinite(out["native_replacement_gc"]) & np.isfinite(out["local_flank_gc"])
        & np.isfinite(out["local_flank_cpg_density"]) & out["local_flank_ambiguous_fraction"].eq(0)
    )
    out["eligible_three_bin_neighborhood"] = False
    for _, group in out.groupby("promoter_id", sort=False):
        ordered = group.sort_values("target_bin_start0", kind="stable")
        starts = ordered["target_bin_start0"].to_numpy(int)
        eligible = ordered["eligible_sequence_target"].to_numpy(bool)
        three = np.zeros(len(ordered), dtype=bool)
        if len(ordered) >= 3:
            consecutive = (starts[1:-1] - starts[:-2] == RESOLUTION) & (starts[2:] - starts[1:-1] == RESOLUTION)
            three[1:-1] = eligible[:-2] & eligible[1:-1] & eligible[2:] & consecutive
        out.loc[ordered.index, "eligible_three_bin_neighborhood"] = three
    return out


def _interval_ambiguous_fraction(fasta, *, start0: int, end0: int) -> float:
    sequence = fasta.fetch(CHROM, int(start0), int(end0)).upper()
    if len(sequence) != int(end0) - int(start0) or not sequence:
        return math.nan
    canonical = sum(b in "ACGT" for b in sequence)
    return float(1.0 - canonical / len(sequence))


def _shared_model_interval(*, tss0: int, site_centers0: list[int], chromosome_length: int) -> tuple[int, int] | None:
    replacement_starts = [int(c) - DONOR_LENGTH // 2 for c in site_centers0]
    replacement_ends = [s + DONOR_LENGTH for s in replacement_starts]
    required_start0 = min([int(tss0) - 2_000, *replacement_starts])
    required_end0 = max([int(tss0) + 2_001, *replacement_ends])
    if required_end0 - required_start0 > MODEL_LENGTH:
        return None
    midpoint = (required_start0 + required_end0) // 2
    start0 = midpoint - MODEL_LENGTH // 2
    start0 = max(0, min(start0, int(chromosome_length) - MODEL_LENGTH))
    end0 = start0 + MODEL_LENGTH
    if not (start0 <= required_start0 and required_end0 <= end0):
        return None
    return start0, end0


def _select_high_low_pairs(contacts: pd.DataFrame, *, fasta, chromosome_length: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    selected_rows: list[pd.Series] = []
    pair_rows: list[dict[str, Any]] = []
    for promoter_id, group in contacts.groupby("promoter_id", sort=False):
        eligible = group[group["eligible_sequence_target"]].copy()
        if eligible.empty:
            continue
        peak_centers = eligible[eligible["eligible_three_bin_neighborhood"]].copy()
        if peak_centers.empty:
            continue
        peak_centers["smooth3_percentile_within_promoter"] = peak_centers["contact_oe_smooth3"].rank(method="average", pct=True)
        peak_centers = peak_centers[
            peak_centers["contact_oe_smooth3"].ge(HIGH_MIN_OE)
            & peak_centers["smooth3_percentile_within_promoter"].ge(HIGH_MIN_PERCENTILE)
        ].sort_values(
            ["contact_oe_smooth3", "contact_oe", "n_nonzero_promoter_bins", "absolute_distance_bp", "target_bin_start0"],
            ascending=[False, False, False, True, True], kind="stable",
        )
        if peak_centers.empty:
            continue
        high: pd.Series | None = None
        primary_lows: list[pd.Series] = []
        seen_high_centers: set[int] = set()
        for _, peak_center in peak_centers.iterrows():
            peak_neighborhood = eligible[
                (eligible["target_bin_center0"].astype(float) - float(peak_center["target_bin_center0"])).abs().le(RESOLUTION)
            ].copy()
            candidate_high = peak_neighborhood.sort_values(
                ["contact_oe", "contact_oe_smooth3", "n_nonzero_promoter_bins", "absolute_distance_bp", "target_bin_start0"],
                ascending=[False, False, False, True, True], kind="stable",
            ).iloc[0].copy()
            high_center0 = int(candidate_high["target_bin_center0"])
            if high_center0 in seen_high_centers:
                continue
            seen_high_centers.add(high_center0)
            if float(candidate_high["contact_oe"]) <= 0:
                continue
            high_flank_gc = float(candidate_high["local_flank_gc"])
            high_flank_cpg = float(candidate_high["local_flank_cpg_density"])
            high_flank_ambiguous = float(candidate_high["local_flank_ambiguous_fraction"])
            if not np.isfinite(high_flank_gc) or not np.isfinite(high_flank_cpg) or high_flank_ambiguous > 0:
                continue
            direction = int(np.sign(float(candidate_high["signed_distance_bp"])))
            distance_difference = (eligible["absolute_distance_bp"].astype(float) - float(candidate_high["absolute_distance_bp"])).abs()
            separated = (eligible["target_bin_center0"].astype(float) - float(candidate_high["target_bin_center0"])).abs()
            local_pool = eligible[
                np.sign(eligible["signed_distance_bp"].astype(float)).eq(direction)
                & distance_difference.le(LOW_DISTANCE_CALIPER) & separated.ge(MIN_SITE_SEPARATION)
            ].copy()
            if len(local_pool) < 2:
                continue
            low_threshold = float(local_pool["contact_oe"].quantile(0.25))
            low_pool = local_pool[
                local_pool["contact_oe"].le(low_threshold) & local_pool["contact_oe"].lt(float(candidate_high["contact_oe"]))
            ].copy()
            low_pool = low_pool[(low_pool["local_flank_gc"].astype(float) - high_flank_gc).abs().le(LOW_GC_TOLERANCE)].copy()
            if len(low_pool) < 2:
                continue
            low_pool["distance_mismatch_scaled"] = (
                low_pool["absolute_distance_bp"].astype(float) - float(candidate_high["absolute_distance_bp"])
            ).abs() / max(float(LOW_DISTANCE_CALIPER), 1.0)
            low_pool["gc_mismatch_scaled"] = (low_pool["local_flank_gc"].astype(float) - high_flank_gc).abs() / max(LOW_GC_TOLERANCE, 1e-6)
            low_pool["cpg_mismatch_scaled"] = (low_pool["local_flank_cpg_density"].astype(float) - high_flank_cpg).abs() / 0.01
            low_pool["match_score"] = low_pool["distance_mismatch_scaled"] + low_pool["gc_mismatch_scaled"] + low_pool["cpg_mismatch_scaled"]
            low_pool = low_pool.sort_values(
                ["match_score", "contact_oe", "absolute_distance_bp", "target_bin_start0"],
                ascending=[True, True, True, True], kind="stable",
            )
            best_primary: tuple[pd.Series, pd.Series] | None = None
            best_primary_key: tuple[float, ...] | None = None
            for combo in combinations([row.copy() for _, row in low_pool.iterrows()], 2):
                centers = [int(row["target_bin_center0"]) for row in combo]
                if abs(centers[0] - centers[1]) < MIN_SITE_SEPARATION:
                    continue
                if _shared_model_interval(
                    tss0=int(candidate_high["tss0"]), site_centers0=[high_center0, *centers],
                    chromosome_length=int(chromosome_length),
                ) is None:
                    continue
                key = (
                    float(sum(float(row["match_score"]) for row in combo)),
                    float(sum(float(row["contact_oe"]) for row in combo)),
                    float(sum(abs(c) for c in centers)),
                )
                if best_primary_key is None or key < best_primary_key:
                    best_primary, best_primary_key = (combo[0].copy(), combo[1].copy()), key
            if best_primary is None:
                continue
            candidate_high["selection_peak_center0"] = int(peak_center["target_bin_center0"])
            candidate_high["selection_peak_oe_mean3"] = float(peak_center["contact_oe_smooth3"])
            candidate_high["selection_peak_percentile"] = float(peak_center["smooth3_percentile_within_promoter"])
            high = candidate_high
            primary_lows = sorted(
                [row.copy() for row in best_primary],
                key=lambda row: (float(row["match_score"]), float(row["contact_oe"]), int(row["target_bin_start0"])),
            )
            break
        if high is None or len(primary_lows) != 2:
            continue
        lows = list(primary_lows)
        high_gc, high_flank_gc = float(high["native_replacement_gc"]), float(high["local_flank_gc"])
        high_flank_cpg, high_flank_ambiguous = float(high["local_flank_cpg_density"]), float(high["local_flank_ambiguous_fraction"])
        chosen_centers = [int(high["target_bin_center0"]), *[int(row["target_bin_center0"]) for row in lows]]
        common_interval = _shared_model_interval(tss0=int(high["tss0"]), site_centers0=chosen_centers, chromosome_length=int(chromosome_length))
        if common_interval is None:
            continue
        model_ambiguous_fraction = _interval_ambiguous_fraction(fasta, start0=int(common_interval[0]), end0=int(common_interval[1]))
        if not np.isfinite(model_ambiguous_fraction):
            continue
        required_coordinates = [
            int(high["tss0"]) - 2_000, int(high["tss0"]) + 2_001,
            *[int(c) - DONOR_LENGTH // 2 for c in chosen_centers],
            *[int(c) - DONOR_LENGTH // 2 + DONOR_LENGTH for c in chosen_centers],
        ]
        model_min_edge_margin = min(
            min(required_coordinates) - int(common_interval[0]), int(common_interval[1]) - max(required_coordinates),
        )
        high["site_role"], high["match_score"] = "high_contact", 0.0
        high["model_interval_start0"], high["model_interval_end0"] = int(common_interval[0]), int(common_interval[1])
        high["model_input_ambiguous_fraction"] = model_ambiguous_fraction
        high["model_min_required_edge_margin_bp"] = int(model_min_edge_margin)
        selected_rows.append(high)
        pair_row: dict[str, Any] = {
            "promoter_id": promoter_id, "gene_id": str(high["gene_id"]), "gene_id_base": str(high["gene_id_base"]),
            "gene_name": str(high["gene_name"]), "gene_type": str(high["gene_type"]), "strand": str(high["strand"]),
            "tss0": int(high["tss0"]), "tss1": int(high["tss1"]),
            "high_target_bin_start0": int(high["target_bin_start0"]), "high_target_center0": int(high["target_bin_center0"]),
            "high_signed_distance_bp": int(high["signed_distance_bp"]), "high_contact_oe": float(high["contact_oe"]),
            "high_contact_oe_smooth3": float(high["contact_oe_smooth3"]),
            "high_selection_peak_center0": int(high["selection_peak_center0"]),
            "high_selection_peak_oe_mean3": float(high["selection_peak_oe_mean3"]),
            "high_selection_peak_percentile": float(high["selection_peak_percentile"]),
            "high_native_replacement_gc": float(high_gc), "high_local_flank_gc": float(high_flank_gc),
            "high_local_flank_cpg_density": float(high_flank_cpg),
            "model_interval_start0": int(common_interval[0]), "model_interval_end0": int(common_interval[1]),
            "model_input_ambiguous_fraction": float(model_ambiguous_fraction),
            "model_min_required_edge_margin_bp": int(model_min_edge_margin),
            "passes_primary_ambiguous_fraction": bool(model_ambiguous_fraction <= 0.01),
            "primary_distance_range_50_500kb": bool(all(50_000 <= int(row["absolute_distance_bp"]) <= 500_000 for row in [high, *primary_lows])),
            "tss_bin_start0": (int(high["tss0"]) // RESOLUTION) * RESOLUTION,
            "contact_edge_cluster": f"{(int(high['tss0']) // RESOLUTION) * RESOLUTION}__{int(high['target_bin_start0'])}",
            "n_low_controls": 2,
        }
        for low_index, low in enumerate(lows, start=1):
            low["site_role"] = f"low_contact_{low_index}"
            low["model_interval_start0"], low["model_interval_end0"] = int(common_interval[0]), int(common_interval[1])
            low["model_input_ambiguous_fraction"] = model_ambiguous_fraction
            low["model_min_required_edge_margin_bp"] = int(model_min_edge_margin)
            selected_rows.append(low)
            prefix = f"low{low_index}"
            pair_row.update({
                f"{prefix}_target_bin_start0": int(low["target_bin_start0"]), f"{prefix}_target_center0": int(low["target_bin_center0"]),
                f"{prefix}_signed_distance_bp": int(low["signed_distance_bp"]), f"{prefix}_contact_oe": float(low["contact_oe"]),
                f"{prefix}_contact_oe_smooth3": float(low["contact_oe_smooth3"]),
                f"{prefix}_native_replacement_gc": float(low["native_replacement_gc"]), f"{prefix}_local_flank_gc": float(low["local_flank_gc"]),
                f"{prefix}_local_flank_cpg_density": float(low["local_flank_cpg_density"]),
                f"{prefix}_absolute_distance_mismatch_bp": abs(int(high["absolute_distance_bp"]) - int(low["absolute_distance_bp"])),
                f"{prefix}_absolute_flank_gc_mismatch": abs(float(high_flank_gc) - float(low["local_flank_gc"])),
                f"{prefix}_absolute_flank_cpg_mismatch": abs(float(high_flank_cpg) - float(low["local_flank_cpg_density"])),
            })
        pair_row["primary_mean_low_contact_oe"] = float(np.mean([float(low["contact_oe"]) for low in primary_lows]))
        pair_row["primary_mean_low_contact_oe_smooth3"] = float(np.mean([float(low["contact_oe_smooth3"]) for low in primary_lows]))
        pair_row["primary_contact_oe_difference"] = float(high["contact_oe"]) - float(pair_row["primary_mean_low_contact_oe"])
        pair_row["primary_contact_oe_smooth3_difference"] = float(high["selection_peak_oe_mean3"]) - float(pair_row["primary_mean_low_contact_oe_smooth3"])
        pair_rows.append(pair_row)
    selected = pd.DataFrame(selected_rows)
    if not selected.empty:
        selected = selected.reset_index(drop=True)
    return selected, pd.DataFrame(pair_rows)


def build_promoter_site_catalog(run_dir: Path, audit: Audit) -> tuple[Path, Path]:
    """Stage 1: chr1 promoter/Hi-C contact catalogue + high/low site
    selection. Zero AlphaGenome cost; deterministic given GENCODE v46 and
    the static 4DNFICSTCJQZ Hi-C file. Checkpointed to run_dir since it is
    ~30-40 minutes of remote Hi-C queries, not an AlphaGenome prediction."""
    out = run_dir / "derived/Figure4EF_promoter_hic_catalog"
    out.mkdir(parents=True, exist_ok=True)
    sites_path, pairs_path = out / "selected_high_low_sites.tsv", out / "selected_high_low_pairs.tsv"
    if sites_path.exists() and pairs_path.exists():
        return sites_path, pairs_path

    import pysam

    ucsc_fasta_path = fetch_ucsc_hg38(run_dir, audit)
    gencode_path = _ensure_gencode(run_dir, audit)
    fasta = pysam.FastaFile(str(ucsc_fasta_path))
    chromosome_length = int(fasta.get_reference_length(CHROM))

    with audit.step("4E/4F: build chr1 GENCODE promoter/exon catalogue"):
        promoters = _load_promoters(gencode_path, chromosome_length=chromosome_length)
        transcript_promoters = _load_transcript_promoters(gencode_path, chromosome_length=chromosome_length)
        merged_starts, merged_ends = _merge_intervals(
            np.concatenate([promoters["promoter_start0"].to_numpy(int), transcript_promoters["promoter_start0"].to_numpy(int)]),
            np.concatenate([promoters["promoter_end0"].to_numpy(int), transcript_promoters["promoter_end0"].to_numpy(int)]),
        )
        merged_exon_starts, merged_exon_ends = _load_merged_exons(gencode_path)
        benchmark = promoters[promoters["benchmark_gene_type"]].copy()

    import hicstraw

    with audit.step(f"4E/4F: query observed HepG2 Hi-C contacts for {len(benchmark)} chr1 protein-coding promoters"):
        hic = hicstraw.HiCFile(HIC_URL)
        zoom = hic.getMatrixZoomData(HIC_CHROM, HIC_CHROM, "oe", "KR", "BP", RESOLUTION)
        chunk_dir = out / "contact_chunks"
        chunk_dir.mkdir(exist_ok=True)
        chunk_size = 25
        frames: list[pd.DataFrame] = []
        for chunk_start in range(0, len(benchmark), chunk_size):
            chunk_number = chunk_start // chunk_size
            chunk_path = chunk_dir / f"contacts_{chunk_number:04d}.parquet"
            if chunk_path.exists():
                frames.append(pd.read_parquet(chunk_path))
                continue
            chunk_frames = [
                _extract_promoter_contacts(
                    zoom, pd.Series(promoter), chromosome_length=chromosome_length,
                    merged_promoter_starts=merged_starts, merged_promoter_ends=merged_ends,
                    merged_exon_starts=merged_exon_starts, merged_exon_ends=merged_exon_ends,
                )
                for promoter in benchmark.iloc[chunk_start: chunk_start + chunk_size].to_dict("records")
            ]
            chunk_frame = pd.concat(chunk_frames, ignore_index=True)
            chunk_frame.to_parquet(chunk_path, index=False)
            frames.append(chunk_frame)
        contacts = pd.concat(frames, ignore_index=True)

    with audit.step("4E/4F: annotate sequence context and select high/low contact-site triples"):
        contacts = _annotate_sequence_context(contacts, fasta=fasta, chromosome_length=chromosome_length)
        selected, pairs = _select_high_low_pairs(contacts, fasta=fasta, chromosome_length=chromosome_length)
        selected.to_csv(sites_path, sep="\t", index=False)
        pairs.to_csv(pairs_path, sep="\t", index=False)
    return sites_path, pairs_path


# --- stage 2: 7-state-per-promoter AlphaGenome transfer scoring -----------

def _load_donors(fasta) -> dict[str, Any]:
    donors = {d.donor_group: d for d in make_asymmetric_315_donors(fasta) if d.donor_group in ("rs127_major", "rs127_minor")}
    return {"major_G": donors["rs127_major"], "minor_T": donors["rs127_minor"]}


def _promoter_order(frame: pd.DataFrame) -> list[str]:
    unique = frame[["promoter_id", "tss0"]].drop_duplicates("promoter_id")
    return unique.sort_values(["tss0", "promoter_id"], kind="stable")["promoter_id"].astype(str).tolist()


def _build_promoter_states(group: pd.DataFrame, *, fasta, donors: dict[str, Any]) -> list[dict[str, Any]]:
    first = group.iloc[0]
    start0, end0 = int(first["model_interval_start0"]), int(first["model_interval_end0"])
    native = fasta.fetch(CHROM, start0, end0).upper()
    if len(native) != MODEL_LENGTH or set(native).difference("ACGTN"):
        raise ValueError(f"Invalid shared FASTA for {first['promoter_id']}")
    tss0 = int(first["tss0"])
    n_fraction = native.count("N") / MODEL_LENGTH
    primary_site_distances = group.loc[group["site_role"].isin(("high_contact", *PRIMARY_LOW_ROLES)), "absolute_distance_bp"].astype(int)
    primary_input_qc = n_fraction <= PRIMARY_MAX_N_FRACTION
    primary_distance = bool((primary_site_distances >= 50_000).all() and (primary_site_distances <= PRIMARY_MAX_DISTANCE).all())
    common = {
        "promoter_id": str(first["promoter_id"]), "gene_id": str(first["gene_id"]), "gene_name": str(first["gene_name"]),
        "gene_type": str(first["gene_type"]), "strand": str(first["strand"]), "tss0": tss0, "tss1": int(first["tss1"]),
        "interval_start0": start0, "interval_end0": end0, "input_n_fraction": n_fraction,
        "primary_input_qc_pass": primary_input_qc, "primary_distance_range": primary_distance,
        "primary_analysis_eligible": primary_input_qc and primary_distance,
    }
    states: list[dict[str, Any]] = [{
        **common, "site_role": "native_shared", "construct": "native", "sequence": native,
        "replacement_start0": None, "replacement_end0": None, "target_bin_center0": None,
    }]
    for role in SITE_ROLES:
        row = group[group["site_role"].astype(str).eq(role)].iloc[0]
        slot_start, slot_end = int(row["replacement_start0"]) - start0, int(row["replacement_end0"]) - start0
        for construct in EDIT_CONSTRUCTS:
            donor = donors[construct]
            sequence = native[:slot_start] + donor.sequence + native[slot_end:]
            if len(sequence) != MODEL_LENGTH:
                raise AssertionError("Replacement changed model-input length")
            states.append({
                **common, "site_role": role, "construct": construct, "sequence": sequence,
                "replacement_start0": int(row["replacement_start0"]), "replacement_end0": int(row["replacement_end0"]),
                "target_bin_center0": int(row["target_bin_center0"]),
            })
    return states


def _select_rna_tracks(metadata: pd.DataFrame, *, strand: str, track_spec: dict[str, Any]) -> np.ndarray:
    """4E/4F only need the RNA_SEQ delta (the Hi-C contact strength -- not
    ATAC/H3K27ac -- is the observed-side variable), so only the strict
    RNA_SEQ selector from the frozen track spec is ported."""
    spec = track_spec["RNA_SEQ"]
    mask = metadata["ontology_curie"].astype(str).eq(HEPG2).to_numpy()
    mask &= metadata["name"].astype(str).isin(spec["names"]).to_numpy()
    mask &= metadata["data_source"].astype(str).eq(spec["data_source"]).to_numpy()
    mask &= metadata["endedness"].astype(str).eq(spec["endedness"]).to_numpy()
    mask &= metadata["strand"].fillna(".").astype(str).isin([strand, "."]).to_numpy()
    indices = np.flatnonzero(mask)
    expected = int(spec["expected_tracks_per_gene_strand"])
    if len(indices) != expected:
        raise RuntimeError(f"Strict RNA_SEQ selector resolved {len(indices)} tracks; expected {expected}")
    return indices


def _summarize_rna(output: Any, *, state: dict[str, Any], track_spec: dict[str, Any]) -> dict[str, Any]:
    genome = _genome_module()
    tss_region = genome.Interval(CHROM, state["tss0"] - TSS_HALF_WIDTH, state["tss0"] + TSS_HALF_WIDTH + 1)
    if output.rna_seq is None:
        raise RuntimeError("Missing requested RNA_SEQ output")
    sliced = output.rna_seq.slice_by_interval(tss_region, match_resolution=True)
    if sliced is None:
        raise RuntimeError("Missing RNA_SEQ slice")
    values = np.asarray(sliced.values, dtype=float)
    if values.ndim == 1:
        values = values[:, np.newaxis]
    metadata = sliced.metadata.copy().reset_index(drop=True)
    indices = _select_rna_tracks(metadata, strand=state["strand"], track_spec=track_spec)
    per_track_mean = np.nanmean(values[:, indices], axis=0)
    return {"rna_hepg2_primary": float(np.nanmean(per_track_mean))}


_GENOME_MODULE_CACHE: list[Any] = []


def _genome_module():
    if not _GENOME_MODULE_CACHE:
        from .figure2 import _ag
        genome, *_ = _ag()
        _GENOME_MODULE_CACHE.append(genome)
    return _GENOME_MODULE_CACHE[0]


def _ag_module():
    from .figure2 import _ag
    return _ag()


def _score_distal_states(run_dir: Path, audit: Audit, states: list[dict[str, Any]], *, batch_size: int, max_workers: int) -> pd.DataFrame:
    """Real, fresh AlphaGenome scoring (RNA_SEQ, HepG2), checkpointed per
    sequence hash -- never seeded from previously-scored data."""
    genome, dna_client, dna_model, _ = _ag_module()
    track_spec = json.loads(TRACK_SPEC_PATH.read_text())
    cache = run_dir / "predictions" / "Figure4EF_state_cache"
    cache.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    missing: list[dict[str, Any]] = []
    for state in states:
        digest = hashlib.sha256(state["sequence"].encode("ascii")).hexdigest()
        state["cache_key"] = digest
        path = cache / digest[:2] / digest[2:4] / f"{digest}.json"
        if path.exists():
            payload = json.loads(path.read_text())
            rows.append({**{k: v for k, v in state.items() if k != "sequence"}, **payload})
        else:
            missing.append(state)
    if missing:
        client = dna_client.create(api_key(), model_version=dna_model.ModelVersion.ALL_FOLDS, timeout=300)
        for start in range(0, len(missing), batch_size):
            batch = missing[start: start + batch_size]
            with audit.step(f"4E/4F: score states {start + 1}-{start + len(batch)} of {len(missing)}"):
                intervals = [genome.Interval(CHROM, s["interval_start0"], s["interval_end0"]) for s in batch]
                outputs = client.predict_sequences(
                    sequences=[s["sequence"] for s in batch],
                    requested_outputs={dna_client.OutputType.RNA_SEQ},
                    ontology_terms=[HEPG2], intervals=intervals, progress_bar=False, max_workers=max_workers,
                )
                for state, output in zip(batch, outputs, strict=True):
                    summary = _summarize_rna(output, state=state, track_spec=track_spec)
                    path = cache / state["cache_key"][:2] / state["cache_key"][2:4] / f"{state['cache_key']}.json"
                    path.parent.mkdir(parents=True, exist_ok=True)
                    path.write_text(json.dumps(summary))
                    rows.append({**{k: v for k, v in state.items() if k != "sequence"}, **summary})
                audit.add_api_calls("4E/4F", len(batch))
                audit.add_api_requests("4E/4F", 1)
    return pd.DataFrame(rows)


def _calculate_effects(predictions: pd.DataFrame, sites: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for promoter_id, site_group in sites.groupby("promoter_id", sort=False):
        pred = predictions[predictions["promoter_id"].astype(str).eq(promoter_id)]
        native_pred = pred[pred["construct"].astype(str).eq("native")]
        if len(native_pred) != 1:
            continue
        native_rna = float(native_pred["rna_hepg2_primary"].iloc[0])
        delta_by_role: dict[str, float] = {}
        for role in SITE_ROLES:
            edited = pred[pred["site_role"].astype(str).eq(role)].set_index("construct")
            if not set(EDIT_CONSTRUCTS).issubset(edited.index):
                delta_by_role = {}
                break
            delta_by_role[role] = float(edited.loc["minor_T", "rna_hepg2_primary"]) - float(edited.loc["major_G", "rna_hepg2_primary"])
        if not delta_by_role:
            continue
        row: dict[str, Any] = {
            "promoter_id": promoter_id, "gene_id": native_pred["gene_id"].iloc[0], "gene_name": native_pred["gene_name"].iloc[0],
            "strand": native_pred["strand"].iloc[0], "tss0": int(native_pred["tss0"].iloc[0]),
            "native_rna": native_rna,
            "primary_input_qc_pass": bool(native_pred["primary_input_qc_pass"].iloc[0]),
            "primary_distance_range": bool(native_pred["primary_distance_range"].iloc[0]),
            "primary_analysis_eligible": bool(native_pred["primary_analysis_eligible"].iloc[0]),
            "high_delta_T_minus_G_rna": delta_by_role["high_contact"],
            "primary_mean_low_delta_T_minus_G_rna": float(np.mean([delta_by_role[r] for r in PRIMARY_LOW_ROLES])),
        }
        row["primary_interaction_rna"] = row["high_delta_T_minus_G_rna"] - row["primary_mean_low_delta_T_minus_G_rna"]
        rows.append(row)
    return pd.DataFrame(rows)


# --- cluster-bootstrap summaries matching Figure 4E/4F source data --------

def _cluster_bootstrap_fraction_positive(
    data: pd.DataFrame, rng: np.random.Generator, n_bootstraps: int,
    value_col: str = "primary_interaction_rna", cluster_col: str = "contact_edge_cluster",
) -> tuple[float, float, float]:
    grouped = data.assign(_positive=data[value_col] > 0).groupby(cluster_col, sort=False)
    cluster_n = grouped.size().to_numpy(dtype=float)
    cluster_positive = grouped["_positive"].sum().to_numpy(dtype=float)
    n_clusters = len(cluster_n)
    point = float(cluster_positive.sum() / cluster_n.sum())
    values = np.empty(n_bootstraps, dtype=float)
    for i in range(n_bootstraps):
        sampled = rng.integers(0, n_clusters, size=n_clusters)
        values[i] = cluster_positive[sampled].sum() / cluster_n[sampled].sum()
    lo, hi = np.quantile(values, [0.025, 0.975])
    return point, float(lo), float(hi)


def _distance_band_table(df: pd.DataFrame, rng: np.random.Generator, n_bootstraps: int) -> pd.DataFrame:
    bins = [50_000, 150_000, 300_000, 500_000, 750_000, 1_000_000]
    labels = ["50–150", "150–300", "300–500", "500–750", "750–1,000"]
    d = df[df["primary_input_qc_pass"].fillna(False)].copy()
    d["absolute_distance_bp"] = d["high_signed_distance_bp"].abs()
    d["distance_label"] = pd.cut(d["absolute_distance_bp"], bins=bins, labels=labels, right=True, include_lowest=True)
    rows = []
    for label in labels:
        subset = d[d["distance_label"] == label].copy()
        point, lo, hi = _cluster_bootstrap_fraction_positive(subset, rng, n_bootstraps)
        rows.append({
            "distance_label": label, "n": len(subset), "n_clusters": subset["contact_edge_cluster"].nunique(),
            "fraction": point, "ci_low": lo, "ci_high": hi,
        })
    return pd.DataFrame(rows)


def _stratified_contact_quintile_table(data: pd.DataFrame, rng: np.random.Generator, n_bootstraps: int) -> pd.DataFrame:
    edges = [50_000, 150_000, 300_000, 500_000, 750_000, 1_000_000]
    labels = ["50–150", "150–300", "300–500", "500–750", "750–1,000"]
    d = data[data["primary_input_qc_pass"].fillna(False)].copy()
    d["absolute_distance_bp"] = d["high_signed_distance_bp"].abs()
    d["distance_stratum"] = pd.cut(d["absolute_distance_bp"], bins=edges, labels=labels, right=True, include_lowest=True)
    rows = []
    for stratum in labels:
        subset = d[d["distance_stratum"] == stratum].copy()
        subset["contact_quintile"] = pd.qcut(subset["primary_contact_oe_difference"], 5, labels=False, duplicates="drop")
        for quintile, qdata in subset.groupby("contact_quintile", sort=True):
            point, lo, hi = _cluster_bootstrap_fraction_positive(qdata, rng, n_bootstraps)
            rows.append({
                "distance_stratum": stratum,
                "distance_midpoint_kb": (edges[labels.index(stratum)] + edges[labels.index(stratum) + 1]) / 2000,
                "contact_quintile": int(quintile) + 1, "n": len(qdata), "n_clusters": qdata["contact_edge_cluster"].nunique(),
                "median_contact_contrast": qdata["primary_contact_oe_difference"].median(),
                "fraction": point, "ci_low": lo, "ci_high": hi,
            })
    return pd.DataFrame(rows)


def run_fig4ef(run_dir: Path, audit: Audit, *, batch_size: int = 128, max_workers: int = 8) -> None:
    import pysam

    sites_path, pairs_path = build_promoter_site_catalog(run_dir, audit)
    ucsc_fasta_path = fetch_ucsc_hg38(run_dir, audit)
    fasta = pysam.FastaFile(str(ucsc_fasta_path))
    sites = pd.read_csv(sites_path, sep="\t")
    pairs = pd.read_csv(pairs_path, sep="\t")

    with audit.step("4E/4F: build 7-state-per-promoter native/T/G transfer design"):
        donors = _load_donors(fasta)
        promoter_ids = _promoter_order(sites)
        all_states: list[dict[str, Any]] = []
        for promoter_id in promoter_ids:
            group = sites[sites["promoter_id"].astype(str).eq(promoter_id)]
            all_states.extend(_build_promoter_states(group, fasta=fasta, donors=donors))

    scored = _score_distal_states(run_dir, audit, all_states, batch_size=batch_size, max_workers=max_workers)
    with audit.step("4E/4F: compute promoter-level T-G high-vs-low interaction and bootstrap summaries"):
        interactions = _calculate_effects(scored, sites)
        merge_keys = ["promoter_id", "gene_id", "gene_name", "strand", "tss0"]
        table = interactions.merge(pairs, on=merge_keys, how="left", validate="one_to_one")
        if table["high_contact_oe"].isna().any():
            raise RuntimeError("4E/4F: contact metadata failed to merge onto scored promoters")
        out = run_dir / "derived/Figure4EF_distal_contact_transfer"
        out.mkdir(parents=True, exist_ok=True)
        table.to_csv(out / "analysis_table.tsv", sep="\t", index=False)

        rng = np.random.default_rng(BOOTSTRAP_SEED)
        bands = _distance_band_table(table, rng, BOOTSTRAP_N)
        bands.to_csv(out / "Figure4E_distance_fraction_positive.tsv", sep="\t", index=False)
        quintiles = _stratified_contact_quintile_table(table, rng, BOOTSTRAP_N)
        quintiles.to_csv(out / "Figure4F_contact_dose_response.tsv", sep="\t", index=False)
        _render4e(bands, run_dir / "figures/Figure4E.svg")
        _render4f(quintiles, run_dir / "figures/Figure4F.svg")


def _render4e(bands: pd.DataFrame, output: Path) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(3.2, 2.4))
    y = np.arange(len(bands))
    x = 100 * bands["fraction"]
    lo = 100 * (bands["fraction"] - bands["ci_low"])
    hi = 100 * (bands["ci_high"] - bands["fraction"])
    ax.errorbar(x, y, xerr=[lo, hi], fmt="o", ms=5, color="#2878B5")
    ax.set_yticks(y, bands["distance_label"])
    ax.invert_yaxis()
    ax.axvline(50, color="#D6D6D6", lw=1.0)
    ax.set_xlabel("Promoters where high-contact\nsite was stronger (%)")
    ax.set_ylabel("Distance to promoter (kb)")
    _save_svg(fig, output)


def _render4f(quintiles: pd.DataFrame, output: Path) -> None:
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(3.6, 2.4))
    norm = mpl.colors.Normalize(vmin=50, vmax=1000)
    cmap = mpl.colormaps["viridis"]
    for stratum, sub in quintiles.groupby("distance_stratum", sort=False):
        color = cmap(norm(float(sub["distance_midpoint_kb"].iloc[0])))
        ax.plot(sub["median_contact_contrast"], 100 * sub["fraction"], "o-", ms=4, color=color, label=f"{stratum} kb")
    ax.axhline(50, color="#D6D6D6", lw=1.0)
    ax.set_xlabel("Contact O/E contrast (high - low)")
    ax.set_ylabel("High-contact stronger (%)")
    ax.legend(frameon=False, fontsize=6, ncol=2)
    _save_svg(fig, output)
