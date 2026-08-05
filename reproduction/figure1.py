"""From-scratch AlphaGenome pipelines for Figure 1B, 1C-middle, and 1E."""

from __future__ import annotations

import gzip
import json
import shutil
import time
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .common import Audit, api_key, download, make_alphagenome_importable, sha256_file, utc_now


VARIANT_CHR = "chr1"
VARIANT_POS = 109_274_968
VARIANT_REF = "G"
VARIANT_ALT = "T"
VARIANT_RSID = "rs12740374"
LIVER = "UBERON:0002107"
HEPG2 = "EFO:0001187"
GENES = ("SORT1", "CELSR2", "PSRC1")
FIG1C_START = 109_209_432
FIG1C_END = 109_340_504

GLGC_URL = "https://csg.sph.umich.edu/willer/public/lipids2013/jointGwasMc_LDL.txt.gz"
CHAIN_URL = "https://imputationserver.sph.umich.edu/resources/chain/hg19_to_hg38.over.chain.gz"
HG38_URL = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/"
    "GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/"
    "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
)
EXPECTED_UNCOMPRESSED = {
    "jointGwasMc_LDL.txt": "387d39882ebb6bd335e7c55e8bce060af3936aa118647278eb4b47f0d67e8e6c",
    "hg38.fa": "9cce8b926416dd96b152deea85188495b75f7ac8d634cc723a017067be8702b7",
}
EXPECTED_CHAIN = "5c0598e500ceb5a78c73086929e8ef993aec309bcafb595139b53d440b125a1d"


def _ag() -> tuple[Any, Any, Any, Any]:
    make_alphagenome_importable()
    from alphagenome.data import genome, ontology
    from alphagenome.models import dna_client, dna_model, variant_scorers

    return genome, ontology, dna_client, (dna_model, variant_scorers)


def _decompress(
    source: Path, destination: Path, expected_sha256: str | None = None
) -> dict[str, object]:
    if not destination.exists():
        partial = destination.with_suffix(destination.suffix + ".part")
        try:
            with gzip.open(source, "rb") as inp, partial.open("wb") as out:
                shutil.copyfileobj(inp, out, length=1024 * 1024)
            partial.replace(destination)
        except Exception:
            partial.unlink(missing_ok=True)
            raise
    observed = sha256_file(destination)
    if expected_sha256 is not None and observed != expected_sha256:
        raise ValueError(
            f"Checksum mismatch for {destination.name}: expected {expected_sha256}, got {observed}"
        )
    return {"path": str(destination), "bytes": destination.stat().st_size, "sha256": observed}


def fetch_fig1c_inputs(run_dir: Path, audit: Audit) -> dict[str, Path]:
    raw = run_dir / "raw"
    glgc_gz = raw / "jointGwasMc_LDL.txt.gz"
    chain = raw / "hg19ToHg38.over.chain.gz"
    fasta_gz = raw / "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
    for url, path in ((GLGC_URL, glgc_gz), (CHAIN_URL, chain), (HG38_URL, fasta_gz)):
        item = download(url, path)
        audit.downloads.append(item)
        audit.save()
    if sha256_file(chain) != EXPECTED_CHAIN:
        raise ValueError("The downloaded UCSC hg19-to-hg38 chain does not match the pinned checksum")
    glgc = raw / "jointGwasMc_LDL.txt"
    fasta = raw / "hg38.fa"
    _decompress(glgc_gz, glgc, EXPECTED_UNCOMPRESSED[glgc.name])
    _decompress(fasta_gz, fasta, EXPECTED_UNCOMPRESSED[fasta.name])
    try:
        import pysam

        if not Path(str(fasta) + ".fai").exists():
            pysam.faidx(str(fasta))
    except ImportError as exc:
        raise RuntimeError("pysam is required to index and validate the downloaded hg38 FASTA") from exc
    return {"glgc": glgc, "chain": chain, "fasta": fasta}


def build_fig1c_variant_table(glgc_path: Path, chain_path: Path, output_path: Path) -> pd.DataFrame:
    """Reconstruct the prespecified 111 variants from the full GLGC file."""
    from pyliftover import LiftOver

    parts: list[pd.DataFrame] = []
    columns = ["SNP_hg19", "rsid", "A1", "A2", "beta", "se", "P-value"]
    for chunk in pd.read_csv(glgc_path, sep="\t", usecols=columns, chunksize=200_000):
        chunk = chunk.rename(columns={"P-value": "p"})
        coordinate = chunk["SNP_hg19"].astype(str).str.extract(r"^chr1:(\d+)$", expand=False)
        keep = coordinate.notna()
        if not keep.any():
            continue
        locus = chunk.loc[keep].copy()
        locus["pos_hg19"] = coordinate.loc[keep].astype(int)
        locus = locus[locus["pos_hg19"].between(109_317_590, 110_317_590)]
        if not locus.empty:
            parts.append(locus)
    if not parts:
        raise ValueError("No chr1 SORT1-locus variants were found in the GLGC input")
    locus = pd.concat(parts, ignore_index=True)
    lift = LiftOver(str(chain_path))
    records: list[dict[str, object]] = []
    for row in locus.itertuples(index=False):
        hits = lift.convert_coordinate("chr1", int(row.pos_hg19) - 1)
        if not hits:
            continue
        chrom, pos0, strand, score = max(hits, key=lambda item: float(item[3]))
        if chrom != "chr1":
            continue
        a1, a2 = str(row.A1).upper(), str(row.A2).upper()
        if strand == "-":
            complement = str.maketrans("ACGT", "TGCA")
            a1, a2 = a1.translate(complement), a2.translate(complement)
        beta = float(row.beta)
        records.append(
            {
                "rsid": str(row.rsid),
                "pos_hg19": int(row.pos_hg19),
                "pos": int(pos0) + 1,
                "effect_allele": a1,
                "other_allele": a2,
                "beta": beta,
                "se": float(row.se),
                "p": float(row.p),
                "liftover_strand": strand,
                "liftover_score": float(score),
                "ldl_raising_allele": a1 if beta >= 0 else a2,
                "ldl_lowering_allele": a2 if beta >= 0 else a1,
            }
        )
    variants = pd.DataFrame(records)
    variants = variants[
        variants["pos"].between(FIG1C_START, FIG1C_END) & variants["p"].lt(5e-8)
    ].copy()
    variants = variants.sort_values(["pos", "rsid"]).drop_duplicates("rsid").reset_index(drop=True)
    variants["pos_mb"] = variants["pos"] / 1e6
    if len(variants) != 111:
        raise ValueError(f"Expected 111 Figure 1C variants after filtering, found {len(variants)}")
    if VARIANT_RSID not in set(variants["rsid"]):
        raise ValueError(f"{VARIANT_RSID} is missing from the reconstructed variant set")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    variants.to_csv(output_path, sep="\t", index=False, float_format="%.10g")
    return variants


def validate_fig1c_reference(table: pd.DataFrame, fasta_path: Path) -> int:
    """Require one reported GLGC allele to match GRCh38 at every lifted site."""
    import pysam

    fasta = pysam.FastaFile(str(fasta_path))
    matched = 0
    for row in table.itertuples(index=False):
        reference = fasta.fetch("chr1", int(row.pos) - 1, int(row.pos)).upper()
        alleles = {str(row.effect_allele).upper(), str(row.other_allele).upper()}
        if reference not in alleles:
            raise ValueError(
                f"Neither GLGC allele matches GRCh38 for {row.rsid} at chr1:{row.pos}"
            )
        matched += 1
    return matched


def _fig1b_ontology_mask(data: Any, term: str) -> np.ndarray:
    if "ontology_curie" not in data.metadata.columns:
        return np.ones(len(data.metadata), dtype=bool)
    return data.metadata.ontology_curie.astype(str).eq(term).to_numpy()


def _fig1b_extract(ref: Any, alt: Any, mask: np.ndarray, label: str) -> tuple[dict[str, np.ndarray], pd.DataFrame]:
    if not mask.any():
        raise RuntimeError(f"No AlphaGenome tracks matched {label}")
    arrays = {
        f"{label}_ref": np.nanmean(np.asarray(ref.values, dtype=np.float32)[:, mask], axis=1),
        f"{label}_alt": np.nanmean(np.asarray(alt.values, dtype=np.float32)[:, mask], axis=1),
        f"{label}_start": np.asarray([int(ref.interval.start)], dtype=np.int64),
        f"{label}_resolution": np.asarray([int(ref.resolution)], dtype=np.int64),
    }
    metadata = ref.metadata.loc[mask].copy()
    metadata.insert(0, "display_track", label)
    metadata.insert(1, "source_track_index", np.flatnonzero(mask))
    return arrays, metadata


def run_fig1b(run_dir: Path, audit: Audit) -> None:
    genome, _, dna_client, model_modules = _ag()
    dna_model, _ = model_modules
    client = dna_client.create(api_key(), model_version=dna_model.ModelVersion.ALL_FOLDS, timeout=300)
    variant = genome.Variant(VARIANT_CHR, VARIANT_POS, VARIANT_REF, VARIANT_ALT, VARIANT_RSID)
    interval = genome.Interval(VARIANT_CHR, VARIANT_POS, VARIANT_POS).resize(
        dna_client.SEQUENCE_LENGTH_500KB
    )
    prediction = client.predict_variant(
        interval=interval,
        variant=variant,
        requested_outputs=[
            dna_client.OutputType.CHIP_TF,
            dna_client.OutputType.DNASE,
            dna_client.OutputType.ATAC,
            dna_client.OutputType.RNA_SEQ,
        ],
        ontology_terms=[LIVER, HEPG2],
    )
    audit.add_api_calls("1B", 1)
    audit.add_api_requests("1B", 1)
    ref = prediction.reference
    alt = prediction.alternate
    factor = ref.chip_tf.metadata.get(
        "transcription_factor", pd.Series(index=ref.chip_tf.metadata.index, dtype=str)
    )
    cebpa = factor.astype(str).str.upper().eq("CEBPA").to_numpy()
    hepatic = _fig1b_ontology_mask(ref.chip_tf, HEPG2)
    if (cebpa & hepatic).any():
        cebpa &= hepatic

    def rna_mask(strand: str) -> np.ndarray:
        tissue = _fig1b_ontology_mask(ref.rna_seq, LIVER)
        if "strand" not in ref.rna_seq.metadata.columns:
            return tissue
        stranded = ref.rna_seq.metadata.strand.astype(str).eq(strand).to_numpy()
        return tissue & stranded if (tissue & stranded).any() else stranded

    specs = [
        ("cebpa", ref.chip_tf, alt.chip_tf, cebpa),
        ("dnase", ref.dnase, alt.dnase, _fig1b_ontology_mask(ref.dnase, LIVER)),
        ("atac", ref.atac, alt.atac, _fig1b_ontology_mask(ref.atac, LIVER)),
        ("rna_plus", ref.rna_seq, alt.rna_seq, rna_mask("+")),
        ("rna_minus", ref.rna_seq, alt.rna_seq, rna_mask("-")),
    ]
    output = run_dir / "predictions" / "Figure1B_locus_tracks"
    output.mkdir(parents=True, exist_ok=True)
    arrays: dict[str, np.ndarray] = {}
    metadata: list[pd.DataFrame] = []
    for label, ref_track, alt_track, mask in specs:
        values, selected = _fig1b_extract(ref_track, alt_track, mask, label)
        arrays.update(values)
        metadata.append(selected)
    np.savez_compressed(output / "tracks.npz", **arrays)
    pd.concat(metadata, ignore_index=True).to_csv(
        output / "selected_track_metadata.tsv", sep="\t", index=False
    )
    (output / "run_metadata.json").write_text(
        json.dumps(
            {
                "created_utc": utc_now(),
                "model_regime": "ALL_FOLDS",
                "input_length_bp": int(dna_client.SEQUENCE_LENGTH_500KB),
                "variant": "GRCh38 chr1:109274968 G>T (rs12740374)",
                "ontology_terms": [LIVER, HEPG2],
                "requested_outputs": ["CHIP_TF", "DNASE", "ATAC", "RNA_SEQ"],
            },
            indent=2,
        )
        + "\n"
    )
    render_fig1b(output / "tracks.npz", run_dir / "figures" / "Figure1B.svg")


def run_fig1e(run_dir: Path, audit: Audit) -> None:
    genome, ontology, dna_client, model_modules = _ag()
    dna_model, _ = model_modules
    client = dna_client.create(api_key(), model_version=dna_model.ModelVersion.FOLD_0, timeout=300)
    variant = genome.Variant(VARIANT_CHR, VARIANT_POS, VARIANT_REF, VARIANT_ALT, VARIANT_RSID)
    interval = genome.Interval(VARIANT_CHR, VARIANT_POS, VARIANT_POS).resize(1_048_576)
    output = client.predict_variant(
        interval=interval,
        variant=variant,
        requested_outputs=[dna_client.OutputType.CONTACT_MAPS],
        ontology_terms=[ontology.from_curie(HEPG2)],
    )
    audit.add_api_calls("1E", 1)
    audit.add_api_requests("1E", 1)
    ref_tracks = output.reference.contact_maps
    alt_tracks = output.alternate.contact_maps
    if ref_tracks is None or alt_tracks is None:
        raise RuntimeError("FOLD_0 CONTACT_MAPS prediction is missing")
    ref = np.nanmean(np.asarray(ref_tracks.values), axis=-1).astype(np.float32)
    alt = np.nanmean(np.asarray(alt_tracks.values), axis=-1).astype(np.float32)
    prediction_dir = run_dir / "predictions" / "Figure1E_fold0_contact"
    prediction_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        prediction_dir / "ref_alt.npz",
        ref=ref,
        alt=alt,
        start=np.asarray(int(ref_tracks.interval.start), dtype=np.int64),
        end=np.asarray(int(ref_tracks.interval.end), dtype=np.int64),
        resolution=np.asarray(int(ref_tracks.resolution), dtype=np.int64),
    )
    centers = int(ref_tracks.interval.start) + (np.arange(ref.shape[0]) + 0.5) * int(
        ref_tracks.resolution
    )
    offsets = (centers - (VARIANT_POS - 1)) / 1000
    indices = np.flatnonzero((offsets >= -60) & (offsets <= 165))
    matrix = ref[np.ix_(indices, indices)]
    positions = offsets[indices]
    table_path = run_dir / "derived" / "Figure1E_fold0_contact.tsv"
    pd.DataFrame(matrix, index=positions, columns=positions).rename_axis(
        index="row_position_from_rs12740374_kb",
        columns="column_position_from_rs12740374_kb",
    ).to_csv(table_path, sep="\t")
    render_fig1e(table_path, run_dir / "figures" / "Figure1E.svg")


def _build_ag_variants(table: pd.DataFrame, fasta_path: Path) -> tuple[list[Any], pd.DataFrame]:
    import pysam

    genome, _, _, _ = _ag()
    fasta = pysam.FastaFile(str(fasta_path))
    variants: list[Any] = []
    rows: list[dict[str, object]] = []
    for row in table.itertuples(index=False):
        ref = fasta.fetch("chr1", int(row.pos) - 1, int(row.pos)).upper()
        alleles = {str(row.effect_allele).upper(), str(row.other_allele).upper()}
        if ref not in alleles:
            raise ValueError(f"Neither GLGC allele matches hg38 for {row.rsid} at chr1:{row.pos}")
        alt = next(allele for allele in alleles if allele != ref)
        variants.append(genome.Variant("chr1", int(row.pos), ref, alt, str(row.rsid)))
        rows.append(
            {
                "rsid": str(row.rsid),
                "pos": int(row.pos),
                "ref_hg38": ref,
                "alt_hg38": alt,
                "ldl_lowering_allele": str(row.ldl_lowering_allele),
                "ldl_raising_allele": str(row.ldl_raising_allele),
            }
        )
    return variants, pd.DataFrame(rows)


def _extract_liver_scores(adata: Any, rsid: str) -> pd.DataFrame:
    gene_rows = adata.obs[adata.obs["gene_name"].isin(GENES)]
    liver_tracks = adata.var[adata.var["ontology_curie"].eq(LIVER)]
    rows: list[dict[str, object]] = []
    for gidx, gene_row in gene_rows.iterrows():
        gene_strand = str(gene_row.get("strand", ""))
        for tidx, track_row in liver_tracks.iterrows():
            track_strand = str(track_row.get("strand", ""))
            if gene_strand in {"+", "-"} and track_strand not in {gene_strand, "."}:
                continue
            rows.append(
                {
                    "rsid": rsid,
                    "gene": str(gene_row["gene_name"]),
                    "gene_id": str(gene_row.get("gene_id", "")),
                    "gene_strand": gene_strand,
                    "ontology_curie": str(track_row["ontology_curie"]),
                    "biosample_name": str(track_row.get("biosample_name", "")),
                    "assay": str(track_row.get("assay", track_row.get("Assay title", ""))),
                    "track_strand": track_strand,
                    "gene_mask_lnfc_ref_to_alt": float(adata[gidx, tidx].X[0, 0]),
                }
            )
    return pd.DataFrame(rows)


def run_fig1c_middle(
    run_dir: Path,
    audit: Audit,
    *,
    batch_size: int,
    max_workers: int,
    max_variants: int | None,
) -> None:
    with audit.step("1C-middle: download GLGC, LiftOver chain, and hg38"):
        inputs = fetch_fig1c_inputs(run_dir, audit)
    variant_path = run_dir / "derived" / "Figure1C_reconstructed_variants.tsv"
    with audit.step("1C-middle: reconstruct the 111-variant set"):
        variants_table = build_fig1c_variant_table(inputs["glgc"], inputs["chain"], variant_path)
    if max_variants is not None:
        causal = variants_table[variants_table["rsid"].eq(VARIANT_RSID)]
        others = variants_table[~variants_table["rsid"].eq(VARIANT_RSID)]
        variants_table = pd.concat([causal, others], ignore_index=True).head(max_variants)
        variants_table = variants_table.sort_values("pos").reset_index(drop=True)
    variants, alleles = _build_ag_variants(variants_table, inputs["fasta"])
    genome, _, dna_client, model_modules = _ag()
    dna_model, variant_scorers = model_modules
    raw_path = run_dir / "predictions" / "Figure1C_middle_raw_track_scores.tsv"
    cached = pd.read_csv(raw_path, sep="\t") if raw_path.exists() else pd.DataFrame()
    complete = set()
    if not cached.empty:
        counts = cached.groupby(["rsid", "gene"]).size().unstack(fill_value=0)
        complete = {
            str(rsid)
            for rsid in counts.index
            if all(gene in counts.columns and counts.loc[rsid, gene] > 0 for gene in GENES)
        }
    variant_by_id = {str(item.name): item for item in variants}
    wanted = [str(value) for value in variants_table["rsid"] if str(value) not in complete]
    client = dna_client.create(api_key(), model_version=dna_model.ModelVersion.ALL_FOLDS, timeout=300)
    interval = genome.Interval("chr1", VARIANT_POS, VARIANT_POS).resize(1_048_576)
    scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"]
    parts = [] if cached.empty else [cached]
    for start in range(0, len(wanted), batch_size):
        batch_ids = wanted[start : start + batch_size]
        with audit.step(f"1C-middle: score variants {start + 1}-{start + len(batch_ids)}"):
            outputs = None
            last_error: Exception | None = None
            for attempt in range(1, 4):
                try:
                    outputs = client.score_variants(
                        intervals=interval,
                        variants=[variant_by_id[value] for value in batch_ids],
                        variant_scorers=[scorer],
                        progress_bar=False,
                        max_workers=max_workers,
                    )
                    break
                except Exception as exc:  # pragma: no cover - live API retry
                    last_error = exc
                    if attempt < 3:
                        time.sleep(2 * attempt)
            if outputs is None:
                raise RuntimeError(f"AlphaGenome scoring batch failed after 3 attempts: {last_error}")
            audit.add_api_calls("1C-middle", len(batch_ids))
            audit.add_api_requests("1C-middle", 1)
            new_parts = []
            for rsid, output_list in zip(batch_ids, outputs, strict=True):
                if not output_list:
                    raise RuntimeError(f"No RNA score object returned for {rsid}")
                extracted = _extract_liver_scores(output_list[0], rsid)
                if set(GENES) - set(extracted["gene"]):
                    raise RuntimeError(f"Missing compatible liver RNA scores for {rsid}")
                new_parts.append(extracted)
            parts.extend(new_parts)
            checkpoint = pd.concat(parts, ignore_index=True).drop_duplicates(
                ["rsid", "gene", "gene_id", "ontology_curie", "biosample_name", "assay", "track_strand"],
                keep="last",
            )
            raw_path.parent.mkdir(parents=True, exist_ok=True)
            checkpoint.to_csv(raw_path, sep="\t", index=False, float_format="%.10g")
            parts = [checkpoint]
    raw = pd.concat(parts, ignore_index=True)
    oriented = raw.merge(alleles, on="rsid", how="left", validate="many_to_one")
    ref_raising = oriented["ref_hg38"].eq(oriented["ldl_raising_allele"])
    alt_lowering = oriented["alt_hg38"].eq(oriented["ldl_lowering_allele"])
    ref_lowering = oriented["ref_hg38"].eq(oriented["ldl_lowering_allele"])
    alt_raising = oriented["alt_hg38"].eq(oriented["ldl_raising_allele"])
    if not ((ref_raising & alt_lowering) | (ref_lowering & alt_raising)).all():
        raise ValueError("Could not orient one or more AlphaGenome scores by LDL allele")
    sign = np.where(ref_raising & alt_lowering, 1.0, -1.0)
    oriented["gene_mask_lnfc_ldl_lowering_minus_raising"] = (
        sign * oriented["gene_mask_lnfc_ref_to_alt"].astype(float)
    )
    aggregate = oriented.groupby(["rsid", "pos", "gene"], as_index=False).agg(
        gene_mask_lnfc=("gene_mask_lnfc_ldl_lowering_minus_raising", "mean"),
        n_compatible_liver_tracks=("gene_mask_lnfc_ldl_lowering_minus_raising", "size"),
    )
    wide = aggregate.pivot(index=["rsid", "pos"], columns="gene", values="gene_mask_lnfc").reset_index()
    wide = wide.rename(columns={gene: f"ag_rna_liver_{gene}" for gene in GENES})
    result = variants_table.merge(wide, on=["rsid", "pos"], how="left", validate="one_to_one")
    result_path = run_dir / "derived" / "Figure1C_middle_ag_scores.tsv"
    result.to_csv(result_path, sep="\t", index=False, float_format="%.10g")
    oriented.to_csv(
        run_dir / "derived" / "Figure1C_middle_oriented_track_scores.tsv",
        sep="\t",
        index=False,
        float_format="%.10g",
    )
    render_fig1c_middle(result_path, run_dir / "figures" / "Figure1C_middle.svg")


def _save_svg(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, bbox_inches="tight", metadata={"Date": "2026-08-02", "Creator": "SORT1-reanalysis"})
    plt.close(fig)


def render_fig1b(source: Path, output: Path) -> None:
    tracks = np.load(source)
    specs = [
        ("cebpa", "CEBPA ChIP"),
        ("dnase", "DNase"),
        ("atac", "ATAC"),
        ("rna_plus", "RNA (+)"),
        ("rna_minus", "RNA (-)"),
    ]
    fig, axes = plt.subplots(5, 1, figsize=(4.1, 5.1), sharex=True)
    for ax, (key, label) in zip(axes, specs, strict=True):
        start = int(tracks[f"{key}_start"][0])
        resolution = int(tracks[f"{key}_resolution"][0])
        positions = start + (np.arange(len(tracks[f"{key}_ref"])) + 0.5) * resolution
        keep = (positions >= 109_250_000) & (positions <= 109_320_000)
        ax.plot(positions[keep], tracks[f"{key}_ref"][keep], color="#1f77b4", lw=0.85)
        ax.plot(positions[keep], tracks[f"{key}_alt"][keep], color="#d62728", lw=0.85)
        ax.axvline(VARIANT_POS, color="goldenrod", lw=0.6)
        ax.set_ylabel(label, fontsize=7)
    axes[-1].set_xlabel("chr1 position (GRCh38)")
    _save_svg(fig, output)


def render_fig1c_middle(source: Path, output: Path) -> None:
    table = pd.read_csv(source, sep="\t")
    fig, axes = plt.subplots(3, 1, figsize=(2.1, 5.0), sharex=True)
    for ax, gene in zip(axes, GENES, strict=True):
        column = f"ag_rna_liver_{gene}"
        causal = table["rsid"].eq(VARIANT_RSID)
        ax.scatter(table.loc[~causal, "pos_mb"], table.loc[~causal, column], s=7, color="#777777")
        ax.scatter(table.loc[causal, "pos_mb"], table.loc[causal, column], s=20, color="#d62728")
        ax.axhline(0, color="black", lw=0.5)
        ax.set_ylabel(f"{gene}\nExon-mask lnFC", fontsize=7)
    axes[0].set_title("AG single variant", fontsize=8)
    axes[-1].set_xlabel("chr1 position (Mb)")
    _save_svg(fig, output)


def render_fig1e(source: Path, output: Path) -> None:
    table = pd.read_csv(source, sep="\t")
    positions = table.iloc[:, 0].to_numpy(float)
    values = table.iloc[:, 1:].to_numpy(float)
    fig, ax = plt.subplots(figsize=(2.5, 2.5))
    image = ax.imshow(
        values,
        origin="lower",
        cmap="inferno",
        vmin=-0.65,
        vmax=2.2,
        extent=(positions[0], positions[-1], positions[0], positions[-1]),
    )
    ax.axvline(0, color="#222222", ls=":", lw=0.75)
    ax.axhline(0, color="#222222", ls=":", lw=0.75)
    ax.set_xlabel("Position from rs12740374 (kb)")
    ax.set_ylabel("Position from rs12740374 (kb)")
    fig.colorbar(image, ax=ax, label="log(O/E)")
    _save_svg(fig, output)
