#!/usr/bin/env python3
"""Write SHA-256 checksums for all currently staged rendered assets."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RENDERED = ROOT / "figures/rendered"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    with (ROOT / "MANIFEST.tsv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    files = []
    for row in rows:
        rel = row.get("release_asset", "")
        if row.get("status") == "exclude" or not rel or rel.startswith("TO_"):
            continue
        path = ROOT / rel
        if path.is_file():
            files.append(path)
        elif path.is_dir():
            files.extend(child for child in path.rglob("*") if child.is_file())
    files = sorted(set(files))
    out = ROOT / "outputs/run_manifests/release_asset_sha256.tsv"
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("path", "size_bytes", "sha256"))
        for path in files:
            writer.writerow((str(path.relative_to(ROOT)), path.stat().st_size, sha256(path)))
    print(f"Recorded {len(files)} release assets in {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
