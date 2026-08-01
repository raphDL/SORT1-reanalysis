#!/usr/bin/env python3
"""Build release-data checksums and an archive-to-release provenance map.

The release tables are allowed to be renamed, but their bytes must match at
least one file in the frozen working archive.  This utility records every
byte-identical archive match so that a release checksum proves lineage rather
than only proving integrity after migration.

Usage:
    python analysis/build_source_provenance.py \
        --workspace-root /path/to/alphaGenome
"""

from __future__ import annotations

import argparse
import csv
import hashlib
from collections import defaultdict
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
SOURCE_DATA = REPO_ROOT / "outputs" / "source_data"
RUN_MANIFESTS = REPO_ROOT / "outputs" / "run_manifests"
GENERATED_IN_RELEASE = {
    "outputs/source_data/Figure1B_locus_tracks/": "analysis/export_fig1b_source.py",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workspace-root",
        type=Path,
        required=True,
        help="Root containing investigation/SORT1_cholesterol_musunuru_2010.",
    )
    args = parser.parse_args()

    archive_root = (
        args.workspace_root.resolve()
        / "investigation"
        / "SORT1_cholesterol_musunuru_2010"
    )
    if not archive_root.is_dir():
        raise SystemExit(f"Frozen archive not found: {archive_root}")

    archive_by_size: dict[int, list[Path]] = defaultdict(list)
    for path in archive_root.rglob("*"):
        if path.is_file():
            archive_by_size[path.stat().st_size].append(path)

    checksum_rows: list[dict[str, object]] = []
    provenance_rows: list[dict[str, object]] = []
    unmatched: list[str] = []

    for release_path in sorted(path for path in SOURCE_DATA.rglob("*") if path.is_file()):
        size = release_path.stat().st_size
        digest = sha256(release_path)
        release_relative = release_path.relative_to(REPO_ROOT).as_posix()
        checksum_rows.append({"path": release_relative, "size_bytes": size, "sha256": digest})

        matches: list[str] = []
        for candidate in archive_by_size[size]:
            if sha256(candidate) == digest:
                matches.append(candidate.relative_to(args.workspace_root.resolve()).as_posix())
        generated_by = next(
            (script for prefix, script in GENERATED_IN_RELEASE.items()
             if release_relative.startswith(prefix)),
            "",
        )
        if not matches and not generated_by:
            unmatched.append(release_relative)
        provenance_rows.append(
            {
                "release_path": release_relative,
                "size_bytes": size,
                "release_sha256": digest,
                "migration_operation": (
                    "generated in release from authorized AlphaGenome API call"
                    if generated_by else "byte-identical copy or rename"
                ),
                "legacy_source_paths": generated_by or ";".join(sorted(matches)),
            }
        )

    if unmatched:
        raise SystemExit("No byte-identical archive source for: " + ", ".join(unmatched))

    RUN_MANIFESTS.mkdir(parents=True, exist_ok=True)
    write_tsv(
        RUN_MANIFESTS / "source_data_sha256.tsv",
        ["path", "size_bytes", "sha256"],
        checksum_rows,
    )
    write_tsv(
        RUN_MANIFESTS / "source_data_provenance.tsv",
        [
            "release_path",
            "size_bytes",
            "release_sha256",
            "migration_operation",
            "legacy_source_paths",
        ],
        provenance_rows,
    )
    print(
        f"Recorded {len(checksum_rows)} release files; archive copies and "
        "explicit release-generated exceptions are fully traced."
    )


if __name__ == "__main__":
    main()
