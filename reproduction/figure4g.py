"""Clean-room AlphaGenome reconstruction for Figure 4G: rs12740374 [G>T]
predicted RNA delta for SORT1/PSRC1/CELSR2 across seven tissue contexts.

Ported from the working archive's
results/figure1/panel_b/heatmap_all_ontologies/run_panel_b_heatmap_all_ontologies.py
(not part of this repository). Unlike 4B-4F, this panel is a single
AlphaGenome variant-scoring call -- one variant, one interval, the RNA_SEQ
RECOMMENDED_VARIANT_SCORERS scorer -- filtered client-side to 11 ontology
tracks and displayed for 7 of them (the archive also scores a second,
unrelated variant, rs646776, used elsewhere and not needed here). See
REPRODUCIBILITY_NEXT_STEPS.md.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .common import Audit, api_key
from .figure2 import CHROM, GENES, RS_POS, _ag

# (display label, ontology_curie) -- the archive's ONTOLOGY_ROWS, minus the
# two rows (hepatocyte, blood) the published panel doesn't display and the
# two specific adipose contexts that only feed the derived "adipose" mean.
ONTOLOGY_ROWS: list[tuple[str, str]] = [
    ("liver", "UBERON:0002107"),
    ("hepatocyte", "CL:0000182"),
    ("adipose subQ", "UBERON:0002190"),
    ("adipose visceral", "UBERON:0010414"),
    ("blood", "UBERON:0000178"),
    ("spleen", "UBERON:0002106"),
    ("lung", "UBERON:0002048"),
    ("brain", "UBERON:0000955"),
    ("heart", "UBERON:0000948"),
    ("kidney", "UBERON:0002113"),
]
ONTOLOGY_CURIE_TO_ROW = {curie: label for label, curie in ONTOLOGY_ROWS}
SELECTED_TERMS = [curie for _, curie in ONTOLOGY_ROWS]
DISPLAY_ROWS = ["liver", "adipose", "spleen", "lung", "brain", "heart", "kidney"]
SEQ_LEN = 2 ** 20


def _extract_scores_by_ontology(adata, genes: tuple[str, ...]) -> pd.DataFrame:
    gene_rows = adata.obs[adata.obs["gene_name"].isin(genes)]
    bio_rows = adata.var[adata.var["ontology_curie"].isin(SELECTED_TERMS)]
    rows: list[dict[str, object]] = []
    for gidx, gene_row in gene_rows.iterrows():
        gene_strand = str(gene_row.get("strand", ""))
        for bidx, bio_row in bio_rows.iterrows():
            track_strand = str(bio_row.get("strand", ""))
            if gene_strand in {"+", "-"} and track_strand not in {gene_strand, "."}:
                continue
            rows.append({
                "gene_name": str(gene_row["gene_name"]),
                "ontology_curie": str(bio_row["ontology_curie"]),
                "raw_score": float(adata[gidx, bidx].X[0, 0]),
            })
    return pd.DataFrame(rows)


def run_fig4g(run_dir: Path, audit: Audit) -> None:
    genome, dna_client, dna_model, variant_scorers = _ag()
    out = run_dir / "derived/Figure4G_tissue_rna"
    out.mkdir(parents=True, exist_ok=True)
    raw_path = run_dir / "predictions/Figure4G_tissue_rna_raw.tsv"
    if raw_path.exists():
        raw = pd.read_csv(raw_path, sep="\t")
    else:
        with audit.step("4G: score rs12740374 RNA_SEQ across tissue ontologies"):
            client = dna_client.create(api_key(), model_version=dna_model.ModelVersion.ALL_FOLDS, timeout=300)
            variant = genome.Variant(chromosome=CHROM, position=RS_POS, reference_bases="G", alternate_bases="T", name="rs12740374")
            interval = genome.Interval(chromosome=CHROM, start=RS_POS, end=RS_POS).resize(SEQ_LEN)
            scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"]
            outputs = client.score_variants(
                intervals=[interval], variants=[variant], variant_scorers=[scorer], progress_bar=False, max_workers=4,
            )
            audit.add_api_calls("4G", 1)
            audit.add_api_requests("4G", 1)
            if not outputs or not outputs[0]:
                raise RuntimeError("No RNA score object returned for rs12740374")
            raw = _extract_scores_by_ontology(outputs[0][0], GENES)
            if raw.empty:
                raise RuntimeError("No target-gene RNA track rows extracted for rs12740374")
            raw_path.parent.mkdir(parents=True, exist_ok=True)
            raw.to_csv(raw_path, sep="\t", index=False)

    with audit.step("4G: aggregate per-context mean RNA delta and derived adipose row"):
        delta = raw.groupby(["gene_name", "ontology_curie"], as_index=False)["raw_score"].mean().rename(columns={"raw_score": "delta"})
        delta["row_label"] = delta["ontology_curie"].map(ONTOLOGY_CURIE_TO_ROW)
        adipose_specific = delta[delta["ontology_curie"].isin(["UBERON:0002190", "UBERON:0010414"])]
        adipose_mean = adipose_specific.groupby("gene_name", as_index=False)["delta"].mean().assign(row_label="adipose")
        delta = pd.concat([delta, adipose_mean], ignore_index=True)
        matrix = delta.pivot_table(index="row_label", columns="gene_name", values="delta", aggfunc="mean").reindex(index=DISPLAY_ROWS, columns=list(GENES))
        table = matrix.reset_index().rename(columns={"row_label": "context"})
        table.to_csv(out / "Figure4G_tissue_rna.tsv", sep="\t", index=False)
