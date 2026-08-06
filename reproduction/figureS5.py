"""Clean-room reproduction for Figure S5 (A-D): held-out FOLD_0 replication
of the SORT1 Kircher MPRA benchmark (Figure 2E), for both substitutions and
single-base deletions.

S5A-C reuse Figure 2E's own ISM-scan machinery (`figure2.py::_score_sort1_
modality`), now parameterized by model version and cache location. The
position range S5 needs (edits measured in both SORT1 and SORT1.2
constructs) is verified to be an exact subset of Figure 2E's own scanned
range (600 positions, identical set), so the ALL_FOLDS pass reuses Figure
2E's own cache directly -- only the FOLD_0 pass is genuinely new. Likewise,
`load_sort1_mpra`'s own `kircher_primary_log2_effect` is verified to equal
the archive's separately-named `kircher_mean_log2_effect` exactly (0 diff)
for every construct measured in both constructs, so it is reused rather
than recomputed from scratch.

S5D is a new pipeline: Kircher single-base deletions on the rs12740374-T
background. ISM does not support deletions, so each edit is scored as an
explicit multi-base Variant spanning [min(deletion, rs12740374),
max(deletion, rs12740374)] that simultaneously encodes the T-background
allele and the deletion, via `client.score_variants` (not ISM). An
rs12740374 G>T baseline (no deletion) is scored once per regime and
subtracted, isolating the deletion-specific effect on top of the
T-background -- both fresh under both regimes, since Figure 2E does not
touch deletions at all.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from .common import Audit, api_key
from .figure2 import (
    BASES, CHROM, GENES, HEPG2, MOTIF_END, MOTIF_START, RS_MAJOR, RS_MINOR, RS_POS,
    _ag, _construct_alleles, _score_sort1_modality, load_sort1_mpra,
)

REGULATORY_WINDOW = 16_384
RNA_WINDOW = 131_072


# --- Shared: the "measured in both forward constructs" edit set ----------

def _both_construct_wide(kircher_path: Path) -> pd.DataFrame:
    """Per-construct SORT1/SORT1.2 raw values, restricted to constructs
    measured in both -- both the edit-set filter and the archive's own
    raw SORT1/SORT1.2 columns come from this same pivot."""
    raw = pd.read_csv(kircher_path, sep="\t", low_memory=False).rename(
        columns={"Chrom": "Chromosome", "Pos": "Position", "Barcodes": "Tags", "Coefficient": "Value", "pValue": "P-Value"}
    )
    raw = raw[
        raw.Release.eq("GRCh38") & raw.Element.isin(("SORT1", "SORT1.2", "SORT1-flip"))
        & raw.Ref.astype(str).str.upper().isin(BASES) & raw.Alt.astype(str).str.upper().isin(BASES)
    ].copy()
    raw.Position = raw.Position.astype(int)
    raw["Value"] = pd.to_numeric(raw["Value"], errors="coerce").round(2)
    raw = _construct_alleles(raw)
    primary = raw[raw.Element.isin(("SORT1", "SORT1.2"))]
    wide = primary.pivot_table(index="variant_id_construct", columns="Element", values="Value", aggfunc="first")
    return wide.dropna(subset=["SORT1", "SORT1.2"])


def _substitution_consensus(kircher_path: Path) -> pd.DataFrame:
    """The 1,790 substitutions measured in both SORT1 and SORT1.2 -- a
    verified subset of Figure 2E's own 1,798-substitution, 600-position
    consensus. `kircher_primary_log2_effect` (Figure 2E's own groupby mean)
    equals the archive's `kircher_mean_log2_effect` exactly for this
    subset (verified: 0 diff), so it is reused, not recomputed."""
    consensus = load_sort1_mpra(kircher_path)
    wide = _both_construct_wide(kircher_path)
    out = consensus[consensus.variant_id_construct.isin(wide.index)].copy()
    out = out.merge(wide[["SORT1", "SORT1.2"]].reset_index(), on="variant_id_construct", how="left")
    out["kircher_mean_log2_effect"] = out["kircher_primary_log2_effect"]
    if len(out) != 1790:
        raise ValueError(f"Expected 1,790 substitutions measured in both constructs, found {len(out)}")
    return out


def run_figs5_substitutions(run_dir: Path, audit: Audit, kircher_path: Path, *, chunk_size: int = 10, max_workers: int = 1) -> None:
    consensus = _substitution_consensus(kircher_path)
    frames: dict[str, pd.DataFrame] = {}
    for regime, cache_dirname, panel in (
        ("ALL_FOLDS", "Figure2E_kircher", "2E"),  # reuses Figure 2E's own cache -- same 600 positions
        ("FOLD_0", "FigureS5_fold0_substitutions", "S5_subs_fold0"),
    ):
        with audit.step(f"S5A-C: score SORT1 Kircher substitutions ({regime})"):
            for modality in ("atac", "h3k27ac", "rna"):
                raw = _score_sort1_modality(
                    consensus, run_dir, audit, modality, chunk_size=chunk_size, max_workers=max_workers,
                    model_version_name=regime, panel=panel, cache_dirname=cache_dirname,
                )
                frames[(regime, modality)] = raw

    keys = ["Position", "GenomeRef", "GenomeAlt", "ConstructRef", "ConstructAlt", "variant_id_genome", "variant_id_construct"]
    combined_parts = []
    for regime in ("ALL_FOLDS", "FOLD_0"):
        matched = consensus.copy()
        for modality in ("atac", "h3k27ac"):
            frame = frames[(regime, modality)]
            extra = [c for c in frame if c not in keys and c not in matched]
            matched = matched.merge(frame[keys + extra].drop_duplicates(keys), on=keys, how="left")
        rna_agg = _aggregate_sort1_rna(frames[(regime, "rna")])
        extra = [c for c in rna_agg if c not in keys and c not in matched]
        matched = matched.merge(rna_agg[keys + extra].drop_duplicates(keys), on=keys, how="left")
        matched["model"] = regime
        matched["in_minor_cebpa_motif"] = matched.Position.between(MOTIF_START, MOTIF_END)
        combined_parts.append(matched)
    combined = pd.concat(combined_parts, ignore_index=True)
    counts = combined.groupby("model")["variant_id_construct"].nunique()
    if counts.nunique() != 1:
        raise ValueError(f"Model comparisons have unequal edit counts: {counts.to_dict()}")

    out_dir = run_dir / "derived"
    out_dir.mkdir(parents=True, exist_ok=True)
    combined.to_csv(out_dir / "FigureS5_substitutions.tsv", sep="\t", index=False)
    _render_model_comparison(combined, "ag_rna_3gene_mean_percent_change", "3-gene mean RNA", run_dir / "figures/FigureS5A.svg")
    _render_model_comparison(combined, "ag_atac_mean_score", "ATAC", run_dir / "figures/FigureS5B.svg")
    _render_model_comparison(combined, "ag_h3k27ac_mean_score", "H3K27ac", run_dir / "figures/FigureS5C.svg")


def _aggregate_sort1_rna(raw: pd.DataFrame) -> pd.DataFrame:
    """Verbatim copy of figure2.py's own aggregation (private there)."""
    keys = ["Position", "GenomeRef", "GenomeAlt", "ConstructRef", "ConstructAlt", "variant_id_genome", "variant_id_construct"]
    per_gene = raw.groupby(keys + ["gene_name"], as_index=False).agg(
        ag_rna_lnfc=("raw_score", "mean"), ag_rna_n_tracks=("raw_score", "size")
    )
    per_gene["ag_rna_percent_change"] = 100.0 * np.expm1(per_gene.ag_rna_lnfc)
    out = per_gene.pivot_table(index=keys, columns="gene_name", values=["ag_rna_lnfc", "ag_rna_n_tracks", "ag_rna_percent_change"]).reset_index()
    out.columns = ["_".join(str(part) for part in col if str(part)).strip("_") if isinstance(col, tuple) else str(col) for col in out.columns]
    for gene in GENES:
        out[f"ag_rna_abs_percent_change_{gene}"] = out[f"ag_rna_percent_change_{gene}"].abs()
    signed = [f"ag_rna_percent_change_{gene}" for gene in GENES]
    absolute = [f"ag_rna_abs_percent_change_{gene}" for gene in GENES]
    out["ag_rna_3gene_mean_percent_change"] = out[signed].mean(axis=1)
    out["ag_rna_3gene_mean_abs_percent_change"] = out[absolute].mean(axis=1)
    out["ag_rna_3gene_max_abs_percent_change"] = out[absolute].max(axis=1)
    return out


def _render_model_comparison(data: pd.DataFrame, column: str, title: str, output: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(6.0, 2.6))
    for ax, model in zip(axes, ("ALL_FOLDS", "FOLD_0"), strict=True):
        sub = data[data.model.eq(model)].dropna(subset=["kircher_mean_log2_effect", column])
        motif = sub.in_minor_cebpa_motif & ~sub.is_rs12740374_position
        rs = sub.is_rs12740374_position
        background = ~(motif | rs)
        ax.scatter(sub.loc[background, "kircher_mean_log2_effect"], sub.loc[background, column], s=4, alpha=0.25, color="#8AA9E8")
        ax.scatter(sub.loc[motif, "kircher_mean_log2_effect"], sub.loc[motif, column], s=6, color="#D94841")
        ax.scatter(sub.loc[rs, "kircher_mean_log2_effect"], sub.loc[rs, column], s=20, marker="D", color="#111111")
        x, y = sub["kircher_mean_log2_effect"].to_numpy(float), sub[column].to_numpy(float)
        r = stats.pearsonr(x, y).statistic
        ax.set_title(f"{model}  r={r:.2f}", fontsize=8)
        ax.set_xlabel("Kircher mean log2 effect", fontsize=7)
    axes[0].set_ylabel(title, fontsize=7)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    plt.close(fig)


# --- S5D: single-base deletions on the rs12740374-T background -----------

def _deletion_consensus(kircher_path: Path) -> pd.DataFrame:
    raw = pd.read_csv(kircher_path, sep="\t", low_memory=False).rename(
        columns={"Chrom": "Chromosome", "Pos": "Position", "Barcodes": "Tags", "Coefficient": "Value", "pValue": "P-Value"}
    )
    raw = raw[
        raw.Release.eq("GRCh38") & raw.Element.isin(("SORT1", "SORT1.2", "SORT1-flip"))
        & raw.Ref.astype(str).str.upper().isin(BASES) & raw.Alt.astype(str).eq("-")
    ].copy()
    raw.Position = raw.Position.astype(int)
    raw["Value"] = pd.to_numeric(raw["Value"], errors="coerce")
    raw["GenomeRef"] = raw.Ref.astype(str).str.upper()
    raw["GenomeAlt"] = ""
    raw["ConstructRef"] = [RS_MINOR if int(pos) == RS_POS else ref for pos, ref in zip(raw.Position, raw.GenomeRef)]
    raw["ConstructAlt"] = ""
    raw["variant_id_genome"] = raw.Position.astype(str) + ":" + raw.GenomeRef + ">DEL"
    raw["variant_id_construct"] = raw.Position.astype(str) + ":" + raw.ConstructRef + ">DEL"
    raw["is_rs12740374_position"] = raw.Position.eq(RS_POS)
    raw["offset"] = raw.Position - RS_POS

    primary = raw[raw.Element.isin(("SORT1", "SORT1.2"))]
    keys = ["Position", "GenomeRef", "GenomeAlt", "ConstructRef", "ConstructAlt", "variant_id_genome", "variant_id_construct", "offset", "is_rs12740374_position"]
    consensus = primary.groupby(keys, as_index=False).agg(
        kircher_primary_log2_effect=("Value", "mean"), kircher_primary_sd=("Value", "std"),
        kircher_primary_n_constructs=("Value", "size"), kircher_primary_min_p=("P-Value", "min"),
        kircher_primary_tags=("Tags", "sum"), kircher_primary_dna=("DNA", "sum"), kircher_primary_rna=("RNA", "sum"),
    ).fillna({"kircher_primary_sd": 0.0}).sort_values("Position").reset_index(drop=True)
    consensus["kircher_primary_abs_log2_effect"] = consensus.kircher_primary_log2_effect.abs()

    wide = primary.pivot_table(index="variant_id_construct", columns="Element", values="Value", aggfunc="first")
    wide = wide.dropna(subset=["SORT1", "SORT1.2"])
    consensus = consensus[consensus.variant_id_construct.isin(set(wide.index))].copy()
    consensus = consensus.merge(wide[["SORT1", "SORT1.2"]].reset_index(), on="variant_id_construct", how="left")
    consensus["kircher_mean_log2_effect"] = consensus["kircher_primary_log2_effect"]
    if len(consensus) != 126:
        raise ValueError(f"Expected 126 single-base deletions measured in both constructs, found {len(consensus)}")
    return consensus


def _build_deletion_variant(row: Any, fasta: Any, genome: Any) -> Any:
    del_pos = int(row.Position)
    if del_pos == RS_POS:
        ref = fasta.fetch(CHROM, del_pos - 1, del_pos).upper()
        return genome.Variant(CHROM, del_pos, ref, "", name=f"{del_pos}:T>DEL_on_T_background")
    start, end = min(del_pos, RS_POS), max(del_pos, RS_POS)
    ref = fasta.fetch(CHROM, start - 1, end).upper()
    if ref[del_pos - start] != str(row.GenomeRef).upper():
        raise ValueError(f"Reference mismatch at {del_pos}: MPRA={row.GenomeRef} FASTA={ref[del_pos - start]}")
    if ref[RS_POS - start] != RS_MAJOR:
        raise ValueError(f"Unexpected rs12740374 reference {ref[RS_POS - start]} in haplotype span")
    alt_list = list(ref)
    alt_list[RS_POS - start] = RS_MINOR
    del alt_list[del_pos - start]
    alt = "".join(alt_list)
    return genome.Variant(CHROM, start, ref, alt, name=f"{del_pos}:{row.ConstructRef}>DEL_on_T_background")


def _score_deletions_regime(consensus: pd.DataFrame, run_dir: Path, audit: Audit, *, model_version_name: str, panel: str, fasta_path: Path) -> pd.DataFrame:
    genome, dna_client, dna_model, variant_scorers = _ag()
    cache_dir = run_dir / "predictions" / f"FigureS5D_{model_version_name.lower()}" / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    import pysam

    fasta = pysam.FastaFile(str(fasta_path))
    positions = sorted(consensus.Position.unique())
    center = (positions[0] + positions[-1]) // 2
    intervals = {
        "regulatory": genome.Interval(CHROM, center, center).resize(REGULATORY_WINDOW),
        "rna": genome.Interval(CHROM, center, center).resize(RNA_WINDOW),
    }
    scorers = {
        "atac": variant_scorers.RECOMMENDED_VARIANT_SCORERS["ATAC"],
        "h3k27ac": variant_scorers.RECOMMENDED_VARIANT_SCORERS["CHIP_HISTONE"],
        "rna": variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"],
    }
    client = dna_client.create(api_key(), model_version=getattr(dna_model.ModelVersion, model_version_name), timeout=300)

    # rs12740374 G>T baseline (no deletion), scored once per regime.
    baseline_path = cache_dir / "baseline.tsv"
    baseline_rna_path = cache_dir / "baseline_rna.tsv"
    if baseline_path.exists() and baseline_rna_path.exists():
        baseline = pd.read_csv(baseline_path, sep="\t").iloc[0]
        baseline_rna = pd.read_csv(baseline_rna_path, sep="\t")
    else:
        with audit.step(f"S5D: score rs12740374 G>T baseline ({model_version_name})"):
            baseline_variant = genome.Variant(CHROM, RS_POS, RS_MAJOR, RS_MINOR, name="rs12740374_G_to_T_baseline")
            atac_out = client.score_variant(interval=intervals["regulatory"], variant=baseline_variant, variant_scorers=[scorers["atac"]])[0]
            h3k27ac_out = client.score_variant(interval=intervals["regulatory"], variant=baseline_variant, variant_scorers=[scorers["h3k27ac"]])[0]
            rna_out = client.score_variant(interval=intervals["rna"], variant=baseline_variant, variant_scorers=[scorers["rna"]])[0]
            baseline_row = {**_chromatin_from_adata(atac_out, "atac"), **_chromatin_from_adata(h3k27ac_out, "h3k27ac")}
            baseline_rna = _rna_rows_from_adata(rna_out)
            pd.DataFrame([baseline_row]).to_csv(baseline_path, sep="\t", index=False)
            baseline_rna.to_csv(baseline_rna_path, sep="\t", index=False)
            audit.add_api_calls(panel, 1)
            audit.add_api_requests(panel, 3)
        baseline = pd.read_csv(baseline_path, sep="\t").iloc[0]

    chunk_size = 24
    chunk_paths = []
    with audit.step(f"S5D: score {len(consensus)} deletions ({model_version_name})"):
        for start_i in range(0, len(consensus), chunk_size):
            chunk = consensus.iloc[start_i : start_i + chunk_size]
            chunk_path = cache_dir / f"chunk_{start_i}_{start_i + len(chunk) - 1}.tsv"
            chunk_paths.append(chunk_path)
            if chunk_path.exists():
                continue
            variants = [_build_deletion_variant(row, fasta, genome) for row in chunk.itertuples(index=False)]
            atac_outputs = client.score_variants(intervals=intervals["regulatory"], variants=variants, variant_scorers=[scorers["atac"]], progress_bar=False, max_workers=1)
            h3k27ac_outputs = client.score_variants(intervals=intervals["regulatory"], variants=variants, variant_scorers=[scorers["h3k27ac"]], progress_bar=False, max_workers=1)
            rna_outputs = client.score_variants(intervals=intervals["rna"], variants=variants, variant_scorers=[scorers["rna"]], progress_bar=False, max_workers=1)
            rows = []
            rna_rows = []
            for row, atac_out, h3k27ac_out, rna_out in zip(chunk.itertuples(index=False), atac_outputs, h3k27ac_outputs, rna_outputs, strict=True):
                record = {"variant_id_construct": row.variant_id_construct}
                record.update(_chromatin_from_adata(atac_out[0], "atac"))
                record.update(_chromatin_from_adata(h3k27ac_out[0], "h3k27ac"))
                rows.append(record)
                rna_frame = _rna_rows_from_adata(rna_out[0])
                rna_frame["variant_id_construct"] = row.variant_id_construct
                rna_rows.append(rna_frame)
            pd.DataFrame(rows).to_csv(chunk_path, sep="\t", index=False)
            pd.concat(rna_rows, ignore_index=True).to_csv(cache_dir / f"chunk_{start_i}_{start_i + len(chunk) - 1}_rna.tsv", sep="\t", index=False)
            audit.add_api_calls(panel, len(chunk))
            audit.add_api_requests(panel, 3)

    scores = pd.concat([pd.read_csv(p, sep="\t") for p in chunk_paths], ignore_index=True)
    rna_raw = pd.concat([pd.read_csv(cache_dir / f"{p.stem}_rna.tsv", sep="\t") for p in chunk_paths], ignore_index=True)

    for track in ("atac", "h3k27ac"):
        mean_col = f"ag_{track}_mean_score"
        scores[f"{mean_col}_combined_vs_hg38"] = scores[mean_col]
        scores[f"ag_{track}_baseline_rs127_g_to_t"] = float(baseline[mean_col])
        scores[mean_col] = scores[f"{mean_col}_combined_vs_hg38"] - float(baseline[mean_col])
        scores[f"ag_{track}_abs_score"] = scores[mean_col].abs()
        scores[f"ag_{track}_max_abs_score"] = scores[mean_col].abs()

    rna_raw["raw_score_combined_vs_hg38"] = rna_raw["raw_score"]
    key_cols = ["gene_name", "gene_id", "assay", "track_strand", "biosample_name"]
    base_rna = baseline_rna[key_cols + ["raw_score"]].rename(columns={"raw_score": "baseline_rs127_g_to_t_raw_score"})
    rna_raw = rna_raw.merge(base_rna, on=key_cols, how="left")
    rna_raw["raw_score"] = rna_raw["raw_score_combined_vs_hg38"] - rna_raw["baseline_rs127_g_to_t_raw_score"]

    merged = consensus.merge(scores, on="variant_id_construct", how="left")
    keys = ["Position", "GenomeRef", "GenomeAlt", "ConstructRef", "ConstructAlt", "variant_id_genome", "variant_id_construct"]
    rna_raw = rna_raw.merge(consensus[keys], on="variant_id_construct", how="left")
    rna_agg = _aggregate_sort1_rna(rna_raw)
    extra = [c for c in rna_agg if c not in keys and c not in merged]
    merged = merged.merge(rna_agg[keys + extra].drop_duplicates(keys), on=keys, how="left")
    merged["model"] = model_version_name
    merged["in_minor_cebpa_motif"] = merged.Position.between(MOTIF_START, MOTIF_END)
    return merged


def _chromatin_from_adata(adata: Any, track: str) -> dict[str, object]:
    meta = adata.var.copy()
    values = np.asarray(adata.X, dtype=float)
    if values.ndim == 1:
        values = values.reshape(1, -1)
    mask = meta.ontology_curie.astype(str).eq(HEPG2).to_numpy()
    if track == "h3k27ac":
        mark = next((c for c in ("histone_mark_code", "histone_mark") if c in meta), None)
        mask = mask & (meta[mark].astype(str).eq("H3K27ac").to_numpy() if mark is not None else np.zeros(len(meta), dtype=bool))
    selected = values[0, mask] if np.any(mask) else np.asarray([], dtype=float)
    prefix = f"ag_{track}"
    return {
        f"{prefix}_n_tracks": int(selected.size),
        f"{prefix}_mean_score": float(np.nanmean(selected)) if selected.size else np.nan,
    }


def _rna_rows_from_adata(adata: Any) -> pd.DataFrame:
    genes = adata.obs[adata.obs.gene_name.isin(GENES)]
    tracks = adata.var[adata.var.ontology_curie.astype(str).eq(HEPG2)]
    rows = []
    for gene_index, gene in genes.iterrows():
        strand = str(gene.get("strand", ""))
        for track_index, track in tracks.iterrows():
            track_strand = str(track.get("strand", ""))
            if strand in {"+", "-"} and track_strand not in {strand, "."}:
                continue
            rows.append({
                "gene_name": str(gene.gene_name), "gene_id": str(gene.get("gene_id", "")),
                "assay": str(track.get("Assay title", track.get("assay", ""))),
                "track_strand": track_strand, "biosample_name": str(track.get("biosample_name", "")),
                "raw_score": float(adata[gene_index, track_index].X[0, 0]),
            })
    return pd.DataFrame(rows)


def run_figs5d(run_dir: Path, audit: Audit, kircher_path: Path, fasta_path: Path) -> None:
    consensus = _deletion_consensus(kircher_path)
    parts = [
        _score_deletions_regime(consensus, run_dir, audit, model_version_name=regime, panel=f"S5D_{regime.lower()}", fasta_path=fasta_path)
        for regime in ("ALL_FOLDS", "FOLD_0")
    ]
    combined = pd.concat(parts, ignore_index=True)
    out = run_dir / "derived/FigureS5D_RNA_deletions.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(out, sep="\t", index=False)
    _render_model_comparison(combined, "ag_rna_3gene_mean_percent_change", "3-gene mean RNA (deletions)", run_dir / "figures/FigureS5D.svg")
