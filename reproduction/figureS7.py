"""Clean-room reproduction for Figure S7: native-locus hotspot audit under
ALL_FOLDS vs FOLD_0.

Reuses Figure 3B's own derived native-locus 501bp ISM outputs
(`native_locus_501bp_all_gene_scores.tsv`, `..._SORT1_hotspots.tsv`) and
Figure 3C's hotspot-window helpers (`_fig3c_windows`,
`_fig3c_max_loss_choices`) to build the same 7 "most-disruptive single-base
edit per position, within an ISM-defined hotspot window" constructs used
throughout Figure 3 (the 6 ISM-ranked hotspots H1-H6 plus the C/EBP control
window), then scores each construct's full 519,488bp T-background sequence
for 3-gene liver RNA(TSS) signal under both ALL_FOLDS and FOLD_0. 8
constructs (native + 7 edited) x 2 models = 16 real AlphaGenome calls,
testing whether held-out-fold predictions agree in sign and rank with the
default model for the manuscript's own mechanistic hotspot claims.
Requires Figure 3B's derived outputs (run panel `3B` first in the same run
directory).
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from .common import Audit, api_key
from .figure2 import CHROM, GENES, LIVER, RS_POS, SEQ_LEN, _ag
from .figure3 import _fig3c_max_loss_choices, _fig3c_windows, _native_t, _save_svg, _summarize

HOTSPOT_ORDER = ("H1", "H2", "H3", "H4", "H5", "H6", "C/EBP")


def _hotspot_designs(native: str, interval: Any, hotspots: pd.DataFrame, all_gene_scores: pd.DataFrame) -> list[dict[str, object]]:
    choices = _fig3c_max_loss_choices(all_gene_scores)
    windows = _fig3c_windows(hotspots)
    designs: list[dict[str, object]] = [
        {"construct": "intact_T", "start_offset": np.nan, "end_offset": np.nan, "sequence": native}
    ]
    for row in windows.itertuples(index=False):
        chars = list(native)
        for offset in range(int(row.start_offset), int(row.end_offset) + 1):
            sequence_index = RS_POS + offset - 1 - int(interval.start)
            chars[sequence_index] = str(choices.loc[offset, "most_disruptive_alt"])
        designs.append({
            "construct": str(row.hotspot), "start_offset": int(row.start_offset), "end_offset": int(row.end_offset),
            "sequence": "".join(chars),
        })
    return designs


def run_figs7(run_dir: Path, audit: Audit, fasta_path: Path) -> None:
    genome, dna_client, dna_model, _ = _ag()
    scores_path = run_dir / "derived/Figure3B_native_501bp_ism/native_locus_501bp_all_gene_scores.tsv"
    hotspots_path = run_dir / "derived/Figure3B_native_501bp_ism/native_locus_501bp_SORT1_hotspots.tsv"
    if not scores_path.exists() or not hotspots_path.exists():
        raise RuntimeError("Figure S7 requires Figure 3B's derived outputs; run panel 3B first (same --run-dir).")

    interval = genome.Interval(CHROM, RS_POS, RS_POS).resize(SEQ_LEN)
    native = _native_t(fasta_path, interval)
    all_gene_scores = pd.read_csv(scores_path, sep="\t")
    hotspots = pd.read_csv(hotspots_path, sep="\t")
    designs = _hotspot_designs(native, interval, hotspots, all_gene_scores)

    cache = run_dir / "predictions/FigureS7_native_hotspot_audit"
    cache.mkdir(parents=True, exist_ok=True)
    frames = []
    for model_name in ("ALL_FOLDS", "FOLD_0"):
        client = None
        rows = []
        with audit.step(f"S7: score native-locus hotspot constructs ({model_name})"):
            for design in designs:
                digest = hashlib.sha256(str(design["sequence"]).encode()).hexdigest()
                path = cache / f"{model_name}_{digest}.tsv"
                if not path.exists():
                    if client is None:
                        client = dna_client.create(api_key(), model_version=getattr(dna_model.ModelVersion, model_name), timeout=300)
                    output = client.predict_sequence(
                        sequence=design["sequence"], requested_outputs={dna_client.OutputType.RNA_SEQ},
                        ontology_terms=[LIVER], interval=interval,
                    )
                    values = _summarize(output, interval)
                    pd.DataFrame([{"gene": g, "rna_tss": values[g]} for g in GENES]).to_csv(path, sep="\t", index=False)
                    audit.add_api_calls("S7", 1)
                    audit.add_api_requests("S7", 1)
                values = pd.read_csv(path, sep="\t").set_index("gene").rna_tss.to_dict()
                for gene in GENES:
                    rows.append({
                        "model": model_name, "construct": design["construct"], "start_offset": design["start_offset"],
                        "end_offset": design["end_offset"], "gene": gene, "rna_tss": values[gene],
                    })
        frame = pd.DataFrame(rows)
        baseline = frame[frame.construct.eq("intact_T")][["gene", "rna_tss"]].rename(columns={"rna_tss": "intact_T_rna_tss"})
        frame = frame.merge(baseline, on="gene", how="left")
        frame["loss_vs_intact_T"] = frame["intact_T_rna_tss"] - frame["rna_tss"]
        frames.append(frame)

    combined = pd.concat(frames, ignore_index=True)
    out = run_dir / "derived/FigureS7_native_sequence_audit.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(out, sep="\t", index=False)
    _render_figs7(combined, run_dir / "figures/FigureS7.svg")


def _render_figs7(scores: pd.DataFrame, path: Path) -> None:
    edited = scores[~scores.construct.eq("intact_T")].copy()
    wide = edited.pivot_table(index=["construct", "gene"], columns="model", values="loss_vs_intact_T").reset_index()
    colors = {name: ("#D94841" if name == "C/EBP" else plt.cm.viridis(i / 7)) for i, name in enumerate(HOTSPOT_ORDER)}
    fig, axes = plt.subplots(1, len(GENES), figsize=(6.9, 2.2))
    for ax, gene in zip(axes, GENES, strict=True):
        subset = wide[wide.gene.eq(gene)].copy()
        for model in ("ALL_FOLDS", "FOLD_0"):
            scale = max(float(subset[model].abs().max()), 1e-12)
            subset[f"{model}_scaled"] = subset[model] / scale
        rho = stats.spearmanr(subset.ALL_FOLDS, subset.FOLD_0)
        ax.axhline(0, color="#666666", lw=0.65)
        ax.axvline(0, color="#666666", lw=0.65)
        for row in subset.itertuples(index=False):
            ax.scatter(row.ALL_FOLDS_scaled, row.FOLD_0_scaled, s=26, color=colors[str(row.construct)], edgecolors="white", linewidths=0.55)
            ax.annotate(str(row.construct), (row.ALL_FOLDS_scaled, row.FOLD_0_scaled), xytext=(3, 2), textcoords="offset points", fontsize=5.6)
        ax.set_xlim(-1.08, 1.08)
        ax.set_ylim(-1.08, 1.08)
        ax.set_aspect("equal")
        ax.set_xlabel("ALL_FOLDS normalized loss", fontsize=7)
        ax.set_title(gene, fontsize=8, fontweight="bold")
        ax.text(0.04, 0.96, f"ρ={rho.statistic:.2f}", transform=ax.transAxes, ha="left", va="top", fontsize=6.5)
    axes[0].set_ylabel("FOLD_0 normalized loss", fontsize=7)
    fig.tight_layout()
    _save_svg(fig, path)
