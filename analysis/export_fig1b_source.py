#!/usr/bin/env python3
"""Export the numerical AlphaGenome track arrays needed to render Figure 1B.

This is the only main-figure analysis for which the legacy working archive
did not retain a compact numerical output. It performs one ALL_FOLDS
``predict_variant`` call and writes a compressed NPZ plus explicit track
selection metadata. The API credential is read only from the
``ALPHAGENOME_API_KEY`` environment variable; it is never written or printed.

After an authorized run, regenerate the source-data checksum/provenance
manifests before attempting the strict public-release gate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

# In an installed release environment, ``alphagenome`` is provided by the
# environment. During release preparation this repository may instead live
# under the original workspace beside ``alphagenome/src``. Discover that
# checkout without embedding a personal absolute path.
for parent in Path(__file__).resolve().parents:
    sdk_source = parent / "alphagenome" / "src"
    if sdk_source.is_dir():
        sys.path.insert(0, str(sdk_source))
        break

from alphagenome.data import genome
from alphagenome.data import track_data as track_data_lib
from alphagenome.models import dna_client, dna_model

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = REPO_ROOT / "outputs" / "source_data" / "Figure1B_locus_tracks"
VARIANT_CHR = "chr1"
VARIANT_POS = 109_274_968
LIVER = "UBERON:0002107"
HEPG2 = "EFO:0001187"


def _ontology_mask(data: track_data_lib.TrackData, term: str) -> np.ndarray:
    if "ontology_curie" not in data.metadata.columns:
        return np.ones(len(data.metadata), dtype=bool)
    return data.metadata.ontology_curie.astype(str).eq(term).to_numpy()


def _cebpa_mask(data: track_data_lib.TrackData) -> np.ndarray:
    factor = data.metadata.get("transcription_factor", pd.Series(index=data.metadata.index, dtype=str))
    mask = factor.astype(str).str.upper().eq("CEBPA").to_numpy()
    hepatic = _ontology_mask(data, HEPG2)
    return mask & hepatic if (mask & hepatic).any() else mask


def _rna_mask(data: track_data_lib.TrackData, strand: str) -> np.ndarray:
    tissue = _ontology_mask(data, LIVER)
    if "strand" not in data.metadata.columns:
        return tissue
    strand_mask = data.metadata.strand.astype(str).eq(strand).to_numpy()
    return tissue & strand_mask if (tissue & strand_mask).any() else strand_mask


def _extract(ref, alt, mask: np.ndarray, label: str) -> tuple[dict[str, np.ndarray], pd.DataFrame]:
    if not mask.any():
        raise RuntimeError(f"No AlphaGenome tracks matched {label}")
    if ref.resolution != alt.resolution or ref.interval != alt.interval:
        raise RuntimeError(f"REF/ALT geometry differs for {label}")
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


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    key = os.environ.get("ALPHAGENOME_API_KEY", "").strip()
    if not key:
        raise RuntimeError("Set ALPHAGENOME_API_KEY for the authorized Figure 1B export")

    client = dna_client.create(
        key,
        model_version=dna_model.ModelVersion.ALL_FOLDS,
        timeout=300,
    )
    variant = genome.Variant(
        chromosome=VARIANT_CHR,
        position=VARIANT_POS,
        reference_bases="G",
        alternate_bases="T",
        name="rs12740374",
    )
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

    specifications = [
        ("cebpa", prediction.reference.chip_tf, prediction.alternate.chip_tf,
         _cebpa_mask(prediction.reference.chip_tf)),
        ("dnase", prediction.reference.dnase, prediction.alternate.dnase,
         _ontology_mask(prediction.reference.dnase, LIVER)),
        ("atac", prediction.reference.atac, prediction.alternate.atac,
         _ontology_mask(prediction.reference.atac, LIVER)),
        ("rna_plus", prediction.reference.rna_seq, prediction.alternate.rna_seq,
         _rna_mask(prediction.reference.rna_seq, "+")),
        ("rna_minus", prediction.reference.rna_seq, prediction.alternate.rna_seq,
         _rna_mask(prediction.reference.rna_seq, "-")),
    ]
    arrays: dict[str, np.ndarray] = {}
    metadata_parts = []
    for label, ref, alt, mask in specifications:
        extracted, metadata = _extract(ref, alt, mask, label)
        arrays.update(extracted)
        metadata_parts.append(metadata)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    npz_path = args.output_dir / "tracks.npz"
    np.savez_compressed(npz_path, **arrays)
    metadata_path = args.output_dir / "selected_track_metadata.tsv"
    pd.concat(metadata_parts, ignore_index=True).to_csv(metadata_path, sep="\t", index=False)
    run = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "model_regime": "ALL_FOLDS",
        "input_length_bp": int(dna_client.SEQUENCE_LENGTH_500KB),
        "variant": "GRCh38 chr1:109274968 G>T (rs12740374)",
        "ontology_terms": [LIVER, HEPG2],
        "requested_outputs": ["CHIP_TF", "DNASE", "ATAC", "RNA_SEQ"],
        "tracks_sha256": hashlib.sha256(npz_path.read_bytes()).hexdigest(),
        "metadata_sha256": hashlib.sha256(metadata_path.read_bytes()).hexdigest(),
    }
    (args.output_dir / "run_metadata.json").write_text(json.dumps(run, indent=2) + "\n")
    print(f"Exported Figure 1B source arrays to {args.output_dir}")


if __name__ == "__main__":
    main()
