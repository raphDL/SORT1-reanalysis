"""Clean-room reproduction for Figure S10 (A-B): distal contact and tissue
context controls.

S10A (50kb-binned contact-associated T-minus-G RNA effect) is a pure
recombination of Figure 4E/4F's own already-computed
`analysis_table.tsv` -- no new AlphaGenome calls at all.

S10B (ALL_FOLDS vs FOLD_0 tissue-context replication) reuses Figure 4G's
own single-variant, all-ontology RNA_SEQ scoring machinery
(`_extract_scores_by_ontology`, `ONTOLOGY_ROWS`, `DISPLAY_ROWS`) under
FOLD_0; if `4G` was already run in the same run directory its ALL_FOLDS
matrix is reused directly (zero new calls), otherwise it is scored fresh
(1 more call). Either way this panel is at most 2 real AlphaGenome calls
total.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from .common import Audit, api_key
from .figure2 import CHROM, GENES, RS_POS, _ag
from .figure4g import DISPLAY_ROWS, ONTOLOGY_CURIE_TO_ROW, SELECTED_TERMS, SEQ_LEN, _extract_scores_by_ontology

DISTANCE_BIN_EDGES = np.arange(0, 1_000_001, 50_000)


# --- S10A: reuses Figure 4E/4F's own analysis_table.tsv -------------------

def run_figs10a(run_dir: Path, audit: Audit) -> None:
    path = run_dir / "derived/Figure4EF_distal_contact_transfer/analysis_table.tsv"
    if not path.exists():
        raise RuntimeError("Figure S10A requires Figure 4E/4F's derived output; run panels 4E,4F first (same --run-dir).")
    with audit.step("S10A: 50kb-binned contact-associated RNA effect (reuses Figure 4E/4F's own scored sequences)"):
        data = pd.read_csv(path, sep="\t")
        data = data[data.primary_input_qc_pass.fillna(False)].copy()
        data["absolute_distance_kb"] = data.high_signed_distance_bp.abs() / 1000
        edges_kb = DISTANCE_BIN_EDGES / 1000
        data["distance_bin_index"] = pd.cut(data.absolute_distance_kb, bins=edges_kb, labels=False, right=True, include_lowest=True)
        rows = []
        for index, (low, high) in enumerate(zip(edges_kb[:-1], edges_kb[1:])):
            subset = data[data.distance_bin_index.eq(index)]
            values = subset.primary_interaction_rna
            rows.append({
                "distance_bin_kb": f"{int(low)}-{int(high)}", "n": int(len(subset)),
                "median_primary_interaction_rna": float(values.median()) if len(subset) else np.nan,
                "q25_primary_interaction_rna": float(values.quantile(0.25)) if len(subset) else np.nan,
                "q75_primary_interaction_rna": float(values.quantile(0.75)) if len(subset) else np.nan,
                "fraction_positive": float((values > 0).mean()) if len(subset) else np.nan,
            })
        summary = pd.DataFrame(rows)
        out = run_dir / "derived/FigureS10_distal_and_tissue_controls"
        out.mkdir(parents=True, exist_ok=True)
        summary.to_csv(out / "FigureS10A_contact_associated_RNA_by_50kb_distance_source.tsv", sep="\t", index=False)
        _render_figs10a(data, summary, run_dir / "figures/FigureS10A.svg")


def _render_figs10a(data: pd.DataFrame, summary: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.0, 2.4))
    midpoints = (DISTANCE_BIN_EDGES[:-1] + DISTANCE_BIN_EDGES[1:]) / 2 / 1000
    ax.scatter(data.absolute_distance_kb, data.primary_interaction_rna, s=4, alpha=0.15, color="#2878B5")
    ax.plot(midpoints, summary.median_primary_interaction_rna, color="black", lw=1.2)
    ax.axhline(0, color="#777777", lw=0.55)
    ax.set_yscale("symlog", linthresh=1e-5)
    ax.set_xlabel("Replacement distance from promoter (kb)")
    ax.set_ylabel("Contact-associated T-G RNA effect")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)


# --- S10B: reuses Figure 4G's own scoring machinery, under FOLD_0 --------

def _tissue_matrix(run_dir: Path, audit: Audit, *, model_version_name: str, panel: str) -> pd.DataFrame:
    genome, dna_client, dna_model, variant_scorers = _ag()
    raw_path = run_dir / f"predictions/FigureS10B_{model_version_name.lower()}_tissue_rna_raw.tsv"
    if raw_path.exists():
        raw = pd.read_csv(raw_path, sep="\t")
    else:
        with audit.step(f"S10B: score rs12740374 RNA_SEQ across tissue ontologies ({model_version_name})"):
            client = dna_client.create(api_key(), model_version=getattr(dna_model.ModelVersion, model_version_name), timeout=300)
            variant = genome.Variant(chromosome=CHROM, position=RS_POS, reference_bases="G", alternate_bases="T", name="rs12740374")
            interval = genome.Interval(chromosome=CHROM, start=RS_POS, end=RS_POS).resize(SEQ_LEN)
            scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"]
            outputs = client.score_variants(intervals=[interval], variants=[variant], variant_scorers=[scorer], progress_bar=False, max_workers=4)
            audit.add_api_calls(panel, 1)
            audit.add_api_requests(panel, 1)
            if not outputs or not outputs[0]:
                raise RuntimeError("No RNA score object returned for rs12740374")
            raw = _extract_scores_by_ontology(outputs[0][0], GENES)
            raw_path.parent.mkdir(parents=True, exist_ok=True)
            raw.to_csv(raw_path, sep="\t", index=False)
    delta = raw.groupby(["gene_name", "ontology_curie"], as_index=False)["raw_score"].mean().rename(columns={"raw_score": "delta"})
    delta["row_label"] = delta.ontology_curie.map(ONTOLOGY_CURIE_TO_ROW)
    adipose_specific = delta[delta.ontology_curie.isin(["UBERON:0002190", "UBERON:0010414"])]
    adipose_mean = adipose_specific.groupby("gene_name", as_index=False)["delta"].mean().assign(row_label="adipose")
    delta = pd.concat([delta, adipose_mean], ignore_index=True)
    matrix = delta.pivot_table(index="row_label", columns="gene_name", values="delta", aggfunc="mean").reindex(index=DISPLAY_ROWS, columns=list(GENES))
    matrix.index.name = "context"
    return matrix


def run_figs10b(run_dir: Path, audit: Audit) -> None:
    all_folds_path = run_dir / "derived/Figure4G_tissue_rna/Figure4G_tissue_rna.tsv"
    if all_folds_path.exists():
        all_folds = pd.read_csv(all_folds_path, sep="\t").set_index("context")[list(GENES)]
    else:
        all_folds = _tissue_matrix(run_dir, audit, model_version_name="ALL_FOLDS", panel="S10B")

    fold0 = _tissue_matrix(run_dir, audit, model_version_name="FOLD_0", panel="S10B")

    with audit.step("S10B: compare ALL_FOLDS and FOLD_0 tissue-context matrices"):
        out = run_dir / "derived/FigureS10_distal_and_tissue_controls"
        out.mkdir(parents=True, exist_ok=True)
        fold0.reset_index().rename(columns={"row_label": "context"}).to_csv(out / "fold0_tissue_matrix.tsv", sep="\t", index=False)

        long_all = all_folds.reset_index().melt("context", var_name="gene", value_name="ALL_FOLDS")
        long_f0 = fold0.reset_index().melt("context", var_name="gene", value_name="FOLD_0")
        merged = long_all.merge(long_f0, on=["context", "gene"], validate="one_to_one")
        merged.to_csv(out / "all_folds_vs_fold0_cells.tsv", sep="\t", index=False)
        x, y = merged.ALL_FOLDS.to_numpy(float), merged.FOLD_0.to_numpy(float)
        summary = {
            "n_cells": int(len(merged)), "pearson_r": float(stats.pearsonr(x, y).statistic),
            "spearman_rho": float(stats.spearmanr(x, y).statistic),
            "sign_agreement_fraction": float((np.sign(x) == np.sign(y)).mean()),
        }
        (out / "S10B_summary.json").write_text(pd.Series(summary).to_json(indent=2) + "\n")
        _render_figs10b(all_folds, fold0, run_dir / "figures/FigureS10B.svg")


def _render_figs10b(all_folds: pd.DataFrame, fold0: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(5.5, 2.65))
    vmax = max(all_folds.to_numpy(dtype=float).__abs__().max(), fold0.to_numpy(dtype=float).__abs__().max())
    for ax, mat, title in zip(axes, [all_folds, fold0], ["ALL_FOLDS", "FOLD_0"], strict=True):
        im = ax.imshow(mat.to_numpy(dtype=float), cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
        ax.set_title(title, fontsize=9)
        ax.set_xticks(range(len(GENES)), GENES, fontsize=7)
        ax.set_yticks(range(len(DISPLAY_ROWS)), DISPLAY_ROWS, fontsize=7)
    fig.colorbar(im, ax=axes, shrink=0.78, label="AlphaGenome RNA delta")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)
