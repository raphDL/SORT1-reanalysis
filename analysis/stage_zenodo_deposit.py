#!/usr/bin/env python3
"""Stage and checksum the four large derived tables for Zenodo upload."""

from __future__ import annotations

import argparse
import csv
import hashlib
import shutil
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
MANIFEST = REPO_ROOT / "outputs" / "run_manifests" / "zenodo_pending_large_outputs.tsv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()

    with MANIFEST.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise RuntimeError("The Zenodo large-output manifest is empty")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    staged = []
    for row in rows:
        source = args.workspace_root / row["legacy_path"]
        if not source.is_file():
            raise FileNotFoundError(source)
        if source.stat().st_size != int(row["size_bytes"]):
            raise RuntimeError(f"Size mismatch before staging: {source}")
        if sha256(source) != row["sha256"]:
            raise RuntimeError(f"SHA-256 mismatch before staging: {source}")
        destination = args.output_dir / source.name
        if not args.verify_only:
            shutil.copy2(source, destination)
        if not destination.is_file() and not args.verify_only:
            raise RuntimeError(f"Copy failed: {destination}")
        if not args.verify_only and sha256(destination) != row["sha256"]:
            raise RuntimeError(f"SHA-256 mismatch after staging: {destination}")
        staged.append({
            "panel": row["panel"],
            "filename": source.name,
            "size_bytes": row["size_bytes"],
            "sha256": row["sha256"],
        })

    if not args.verify_only:
        manifest = args.output_dir / "SHA256SUMS.tsv"
        with manifest.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=staged[0].keys(), delimiter="\t")
            writer.writeheader()
            writer.writerows(staged)
        print(f"Staged {len(staged)} files in {args.output_dir}")
    else:
        print(f"Verified {len(staged)} source files; nothing copied")


if __name__ == "__main__":
    main()
