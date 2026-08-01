#!/usr/bin/env python3
"""Freeze the author-approved manuscript package into the release repository.

The source directory must contain the manuscript text, two supplementary
tables, four assembled main-figure SVGs and ten assembled supplementary SVGs.
Files are copied byte-for-byte under stable release names and checksummed.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT_DIR = ROOT / "manuscript"
ASSEMBLED_DIR = ROOT / "figures" / "assembled"

SOURCE_TO_RELEASE = {
    "SORT1 Manuscript.txt": MANUSCRIPT_DIR / "SORT1_Manuscript.txt",
    "Supplementary_Table_1.tsv": MANUSCRIPT_DIR / "Supplementary_Table_1.tsv",
    "Supplementary_Table_2.tsv": MANUSCRIPT_DIR / "Supplementary_Table_2.tsv",
    "Figure1.svg": ASSEMBLED_DIR / "Figure1.svg",
    "figure2.svg": ASSEMBLED_DIR / "Figure2.svg",
    "figure3.svg": ASSEMBLED_DIR / "Figure3.svg",
    "figure4.svg": ASSEMBLED_DIR / "Figure4.svg",
    **{
        f"S{index}.svg": ASSEMBLED_DIR / f"FigureS{index}.svg"
        for index in range(1, 11)
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    args = parser.parse_args()
    source_dir = args.source_dir.resolve()

    missing = [name for name in SOURCE_TO_RELEASE if not (source_dir / name).is_file()]
    if missing:
        raise SystemExit("Missing submission assets: " + ", ".join(missing))

    MANUSCRIPT_DIR.mkdir(parents=True, exist_ok=True)
    ASSEMBLED_DIR.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for source_name, destination in SOURCE_TO_RELEASE.items():
        source = source_dir / source_name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
        rows.append(
            {
                "release_path": destination.relative_to(ROOT).as_posix(),
                "source_filename": source_name,
                "size_bytes": destination.stat().st_size,
                "sha256": sha256(destination),
            }
        )

    checksum_path = MANUSCRIPT_DIR / "SUBMISSION_SNAPSHOT_SHA256.tsv"
    with checksum_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("release_path", "source_filename", "size_bytes", "sha256"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)

    metadata = {
        "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_directory_name": source_dir.name,
        "main_figures": [f"Figure{i}" for i in range(1, 5)],
        "supplementary_figures": [f"FigureS{i}" for i in range(1, 11)],
        "supplementary_numbering_final": True,
        "copy_operation": "byte-identical copy with filename normalization",
        "checksum_manifest": checksum_path.relative_to(ROOT).as_posix(),
    }
    (MANUSCRIPT_DIR / "SUBMISSION_SNAPSHOT.json").write_text(
        json.dumps(metadata, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"Frozen {len(rows)} submission assets; numbering is Figure S1-S10.")


if __name__ == "__main__":
    main()

