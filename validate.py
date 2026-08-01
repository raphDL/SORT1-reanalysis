#!/usr/bin/env python3
"""Validate the compact SORT1 reanalysis release.

The default validation checks repository structure and the panel manifest.
During migration, ``--workspace-root`` additionally verifies that every
non-placeholder legacy asset still exists in the frozen working archive.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REQUIRED_PATHS = (
    "README.md",
    "MANIFEST.tsv",
    "MANIFEST_NOTES.md",
    "PRESPECIFIED.md",
    "METHOD_CONVENTIONS.md",
    "LICENSE",
    "CITATION.cff",
    "RELEASE_CHECKLIST.md",
    "AUDIT_REPORT.md",
    "environment.yml",
    "config.yaml",
    "data/SOURCES.tsv",
    "analysis",
    "figures",
    "outputs/source_data",
    "outputs/run_manifests",
    "outputs/run_manifests/zenodo_deposit_status.tsv",
    "outputs/manuscript_results.tsv",
)
ALLOWED_STATUS = {
    "final",
    "provisional-final",
    "provisional",
    "author-layout",
    "exclude",
}
ALLOWED_MODEL_REGIME_TOKENS = {
    "NA",
    "EXPERIMENTAL",
    "ALL_FOLDS",
    "FOLD_0",
    "MATCHED_HELDOUT",
}
PLACEHOLDER_PREFIXES = ("TO_ADD", "TO_IDENTIFY", "TO_FILL", "TO_REGENERATE")
PERSONAL_PATH_MARKERS = ("/Users/", "/home/", "C:\\Users\\")
SCRIPT_DIRS = ("analysis", "figures", "src")
SCRIPT_SUFFIXES = {".py", ".sh", ".zsh", ".R", ".yaml", ".yml", ".toml"}
PUBLIC_SOURCE_REQUIRED_FIELDS = (
    "provider_or_citation",
    "accession_or_url",
    "genome_build",
    "access_date",
    "checksum",
    "license_or_terms",
    "used_in",
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def check_source_data_checksums(root: Path) -> list[str]:
    """Item 9: source-table and manifest hashes must match."""
    errors: list[str] = []
    manifest_path = root / "outputs/run_manifests/source_data_sha256.tsv"
    if not manifest_path.exists():
        return [f"missing outputs/run_manifests/source_data_sha256.tsv"]
    recorded = read_tsv(manifest_path)
    recorded_paths = set()
    for row in recorded:
        rel = row["path"]
        recorded_paths.add(rel)
        target = root / rel
        if not target.exists():
            errors.append(f"source_data_sha256.tsv: missing file: {rel}")
            continue
        if sha256(target) != row["sha256"]:
            errors.append(f"source_data_sha256.tsv: checksum mismatch: {rel}")
    source_data = root / "outputs/source_data"
    if source_data.exists():
        for path in sorted(source_data.rglob("*")):
            if path.is_file():
                rel = str(path.relative_to(root))
                if rel not in recorded_paths:
                    errors.append(f"source_data_sha256.tsv: unrecorded file present: {rel}")
    return errors


def check_source_data_provenance(root: Path) -> list[str]:
    """Every compact release table must trace to an archived byte-identical source."""
    errors: list[str] = []
    path = root / "outputs/run_manifests/source_data_provenance.tsv"
    if not path.exists():
        return ["missing outputs/run_manifests/source_data_provenance.tsv"]
    rows = read_tsv(path)
    mapped = {row.get("release_path", "") for row in rows}
    for row in rows:
        rel = row.get("release_path", "")
        digest = row.get("release_sha256", "")
        target = root / rel
        if not rel or not target.is_file():
            errors.append(f"source_data_provenance.tsv: absent release file: {rel or '<blank>'}")
        elif sha256(target) != digest:
            errors.append(f"source_data_provenance.tsv: checksum mismatch: {rel}")
        if not row.get("legacy_source_paths", "").strip():
            errors.append(f"source_data_provenance.tsv: no archived source recorded: {rel}")
    for target in sorted((root / "outputs/source_data").rglob("*")):
        if target.is_file():
            rel = str(target.relative_to(root))
            if rel not in mapped:
                errors.append(f"source_data_provenance.tsv: unrecorded source-data file: {rel}")
    return errors


def check_analysis_parameters(rows: list[dict[str, str]], root: Path) -> list[str]:
    """All selected AlphaGenome panels must record their model input geometry."""
    errors: list[str] = []
    path = root / "outputs/run_manifests/panel_analysis_parameters.tsv"
    if not path.exists():
        return ["missing outputs/run_manifests/panel_analysis_parameters.tsv"]
    parameters = {row["panel"]: row for row in read_tsv(path)}
    for row in rows:
        if row["status"] in {"exclude", "author-layout"}:
            continue
        tokens = {token.strip() for token in row["model_regime"].split(";")}
        if tokens <= {"NA", "EXPERIMENTAL", ""}:
            continue
        panel = row["panel"]
        if panel not in parameters:
            errors.append(f"panel_analysis_parameters.tsv: missing AlphaGenome panel {panel}")
            continue
        record = parameters[panel]
        for field in (
            "model_regime", "prediction_context_bp", "output_modality",
            "ontology_or_biosample", "track_selection", "scorer_or_measure",
            "aggregation_window_or_endpoint", "output_summary", "parameter_source",
        ):
            if not record.get(field, "").strip():
                errors.append(f"panel_analysis_parameters.tsv: {panel}: blank {field}")
    return errors


def check_results_ledger(root: Path) -> list[str]:
    errors: list[str] = []
    path = root / "outputs/manuscript_results.tsv"
    if not path.exists():
        return ["missing outputs/manuscript_results.tsv"]
    rows = read_tsv(path)
    if not rows:
        return ["outputs/manuscript_results.tsv contains no results"]
    seen: set[str] = set()
    for row in rows:
        result_id = row.get("result_id", "")
        if not result_id:
            errors.append("manuscript_results.tsv: blank result_id")
        elif result_id in seen:
            errors.append(f"manuscript_results.tsv: duplicate result_id: {result_id}")
        seen.add(result_id)
        source = row.get("source_table", "")
        if not source or not (root / source).exists():
            errors.append(f"manuscript_results.tsv: {result_id}: absent source_table: {source or '<blank>'}")
        if row.get("status") != "verified_from_release_table":
            errors.append(f"manuscript_results.tsv: {result_id}: unverified status {row.get('status', '')!r}")
    return errors


def check_no_personal_paths(root: Path) -> list[str]:
    """Item 9: scripts must contain no personal absolute paths."""
    errors: list[str] = []
    for dirname in SCRIPT_DIRS:
        directory = root / dirname
        if not directory.exists():
            continue
        for path in sorted(directory.rglob("*")):
            if not path.is_file() or path.suffix not in SCRIPT_SUFFIXES:
                continue
            text = path.read_text(encoding="utf-8", errors="replace")
            for marker in PERSONAL_PATH_MARKERS:
                if marker in text:
                    errors.append(
                        f"{path.relative_to(root)}: contains personal absolute path marker {marker!r}"
                    )
    return errors


def check_zenodo_pending_hashes(root: Path) -> list[str]:
    """Large outputs excluded from git must still carry a real SHA-256, not a placeholder."""
    errors: list[str] = []
    path = root / "outputs/run_manifests/zenodo_pending_large_outputs.tsv"
    if not path.exists():
        return errors
    for row in read_tsv(path):
        digest = row.get("sha256", "")
        if not digest or digest.startswith(PLACEHOLDER_PREFIXES) or len(digest) != 64:
            errors.append(f"zenodo_pending_large_outputs.tsv: invalid sha256 for {row.get('panel', '?')}: {digest!r}")
    return errors


def check_zenodo_deposit(root: Path, strict: bool) -> tuple[list[str], list[str]]:
    """Require a real deposited record only at the strict public-release gate."""
    errors: list[str] = []
    warnings: list[str] = []
    pending = root / "outputs/run_manifests/zenodo_pending_large_outputs.tsv"
    status_path = root / "outputs/run_manifests/zenodo_deposit_status.tsv"
    if not status_path.exists():
        return (["missing outputs/run_manifests/zenodo_deposit_status.tsv"] if strict else [], [])
    rows = read_tsv(status_path)
    if len(rows) != 1:
        return (["zenodo_deposit_status.tsv must contain exactly one status row"], [])
    row = rows[0]
    recorded_manifest_hash = row.get("large_output_manifest_sha256", "")
    if pending.exists() and recorded_manifest_hash != sha256(pending):
        errors.append("zenodo_deposit_status.tsv: large-output manifest checksum mismatch")
    try:
        manifest_rows = read_tsv(pending)
        expected_count = len(manifest_rows)
        expected_size = sum(int(item["size_bytes"]) for item in manifest_rows)
        if int(row.get("file_count", "-1")) != expected_count:
            errors.append("zenodo_deposit_status.tsv: file_count does not match manifest")
        if int(row.get("total_size_bytes", "-1")) != expected_size:
            errors.append("zenodo_deposit_status.tsv: total_size_bytes does not match manifest")
    except (ValueError, FileNotFoundError):
        errors.append("zenodo_deposit_status.tsv: invalid count/size metadata")
    deposited = row.get("status", "").lower() == "deposited"
    doi = row.get("doi", "")
    record_url = row.get("record_url", "")
    valid_doi = bool(re.fullmatch(r"10\.5281/zenodo\.\d+", doi))
    valid_url = bool(re.fullmatch(r"https://zenodo\.org/records/\d+", record_url))
    if not (deposited and valid_doi and valid_url):
        message = "Zenodo deposit is pending or lacks a valid DOI/record URL"
        (errors if strict else warnings).append(message)
    return errors, warnings


def check_release_asset_checksums(root: Path) -> list[str]:
    errors: list[str] = []
    path = root / "outputs/run_manifests/release_asset_sha256.tsv"
    if not path.exists():
        return ["missing outputs/run_manifests/release_asset_sha256.tsv"]
    recorded: set[str] = set()
    for row in read_tsv(path):
        rel = row.get("path", "")
        recorded.add(rel)
        target = root / rel
        if not target.is_file():
            errors.append(f"release_asset_sha256.tsv: absent asset: {rel}")
            continue
        if str(target.stat().st_size) != row.get("size_bytes", ""):
            errors.append(f"release_asset_sha256.tsv: size mismatch: {rel}")
        if sha256(target) != row.get("sha256", ""):
            errors.append(f"release_asset_sha256.tsv: checksum mismatch: {rel}")
        prefix = target.read_bytes()[:8]
        if target.suffix.lower() == ".pdf" and not prefix.startswith(b"%PDF-"):
            errors.append(f"release asset extension/content mismatch: {rel} is not PDF")
        if target.suffix.lower() == ".svg" and b"<svg" not in target.read_bytes()[:1024]:
            errors.append(f"release asset extension/content mismatch: {rel} is not SVG")
    # PNG/PDF convenience exports made by deterministic renderers may be
    # git-ignored. Only assets selected by MANIFEST.tsv are release artifacts.
    selected_assets = {
        row.get("release_asset", "")
        for row in read_tsv(root / "MANIFEST.tsv")
        if row.get("status") != "exclude" and row.get("release_asset", "")
    }
    for rel in sorted(selected_assets):
        target = root / rel
        if target.is_file() and rel not in recorded:
            errors.append(f"release_asset_sha256.tsv: unrecorded selected asset: {rel}")
        elif target.is_dir():
            for child in sorted(target.rglob("*")):
                if child.is_file():
                    child_rel = str(child.relative_to(root))
                    if child_rel not in recorded:
                        errors.append(f"release_asset_sha256.tsv: unrecorded component asset: {child_rel}")
    return errors


def check_rows(rows: list[dict[str, str]], root: Path, workspace_root: Path | None,
                require_release_assets: bool) -> tuple[list[str], list[str]]:
    errors: list[str] = []
    warnings: list[str] = []

    panels = [row["panel"] for row in rows]
    duplicate_panels = sorted(panel for panel, n in Counter(panels).items() if n > 1)
    if duplicate_panels:
        errors.append(f"duplicate panel identifiers: {', '.join(duplicate_panels)}")

    zenodo_manifest = root / "outputs/run_manifests/zenodo_pending_large_outputs.tsv"
    zenodo_panels = {
        item["panel"] for item in read_tsv(zenodo_manifest)
    } if zenodo_manifest.exists() else set()

    for row in rows:
        panel = row["panel"]
        status = row["status"]
        if status not in ALLOWED_STATUS:
            errors.append(f"{panel}: invalid status {status!r}")

        selected = status != "exclude"
        if selected and not row["analysis_id"]:
            errors.append(f"{panel}: missing analysis_id")

        if status not in {"author-layout", "exclude"} and not row["release_source_table"]:
            errors.append(f"{panel}: missing planned release source table")

        if require_release_assets and selected:
            required_columns = ["release_asset"]
            if status != "author-layout":
                required_columns.append("release_source_table")
            for column in required_columns:
                value = row[column]
                if not value or value.startswith(PLACEHOLDER_PREFIXES):
                    errors.append(f"{panel}: unresolved {column}: {value or '<blank>'}")
                elif not (root / value).exists():
                    if column == "release_source_table" and panel in zenodo_panels:
                        continue
                    errors.append(f"{panel}: absent {column}: {value}")

        legacy = row["legacy_asset"]
        if workspace_root and legacy and not legacy.startswith(PLACEHOLDER_PREFIXES):
            if not (workspace_root / legacy).exists():
                errors.append(f"{panel}: legacy asset not found: {legacy}")

        if "TO_" in "\t".join(row.values()) and selected:
            warnings.append(f"{panel}: contains a migration placeholder")

        # Item 9: final programmatic panels should have a rendering command.
        if status in {"final", "provisional-final"} and status != "author-layout":
            asset = row["release_asset"]
            if not asset or asset.startswith(PLACEHOLDER_PREFIXES) or not (root / asset).exists():
                warnings.append(f"{panel}: final/provisional-final panel has no rendered release_asset yet")

        # Item 9: manual panels must list their component assets / adaptation source.
        if status == "author-layout":
            if not row["notes"].strip():
                errors.append(f"{panel}: author-layout panel has no notes describing composition")
            if row["external_sources"].strip() and not row["notes"].strip():
                errors.append(f"{panel}: adapted panel missing citation/adaptation notes")

        # Item 9: AlphaGenome-derived analyses must declare a model regime.
        if row["model_regime"]:
            tokens = {t.strip() for t in row["model_regime"].split(";") if t.strip()}
            unknown = tokens - ALLOWED_MODEL_REGIME_TOKENS
            if unknown:
                errors.append(f"{panel}: unrecognized model_regime token(s): {', '.join(sorted(unknown))}")

    return errors, warnings


def check_sources_tsv(root: Path, strict: bool) -> tuple[list[str], list[str]]:
    """Public datasets must have complete, explicit provenance metadata."""
    errors: list[str] = []
    warnings: list[str] = []
    path = root / "data/SOURCES.tsv"
    if not path.exists():
        return errors, warnings
    for row in read_tsv(path):
        for field in PUBLIC_SOURCE_REQUIRED_FIELDS:
            value = row.get(field, "").strip()
            if not value or value.startswith(PLACEHOLDER_PREFIXES):
                message = f"data/SOURCES.tsv: {row.get('source_id', '?')}: unresolved {field} ({value or '<blank>'})"
                (errors if strict else warnings).append(message)
    return errors, warnings


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--workspace-root",
        type=Path,
        help="Optional alphaGenome workspace root used to check legacy paths.",
    )
    parser.add_argument(
        "--require-release-assets",
        action="store_true",
        help="Fail when selected release assets or source tables are absent.",
    )
    args = parser.parse_args()

    errors: list[str] = []
    warnings: list[str] = []

    for relative in REQUIRED_PATHS:
        if not (ROOT / relative).exists():
            errors.append(f"missing required path: {relative}")

    manifest_path = ROOT / "MANIFEST.tsv"
    if not manifest_path.exists():
        for message in errors:
            print(f"ERROR: {message}")
        return 1

    rows = read_tsv(manifest_path)

    checksum_path = ROOT / "outputs/run_manifests/SHA256SUMS"
    if checksum_path.exists():
        for line_number, line in enumerate(
            checksum_path.read_text(encoding="utf-8").splitlines(), start=1
        ):
            if not line.strip():
                continue
            try:
                expected, filename = line.split(maxsplit=1)
            except ValueError:
                errors.append(f"SHA256SUMS:{line_number}: malformed line")
                continue
            target = checksum_path.parent / filename.strip()
            if not target.exists():
                errors.append(f"SHA256SUMS: missing file: {filename.strip()}")
            elif sha256(target) != expected:
                errors.append(f"SHA256SUMS: checksum mismatch: {filename.strip()}")

    row_errors, row_warnings = check_rows(
        rows, ROOT, args.workspace_root, args.require_release_assets
    )
    errors.extend(row_errors)
    warnings.extend(row_warnings)

    errors.extend(check_source_data_checksums(ROOT))
    errors.extend(check_source_data_provenance(ROOT))
    errors.extend(check_analysis_parameters(rows, ROOT))
    errors.extend(check_results_ledger(ROOT))
    errors.extend(check_no_personal_paths(ROOT))
    errors.extend(check_zenodo_pending_hashes(ROOT))
    zenodo_errors, zenodo_warnings = check_zenodo_deposit(ROOT, args.require_release_assets)
    errors.extend(zenodo_errors)
    warnings.extend(zenodo_warnings)
    errors.extend(check_release_asset_checksums(ROOT))
    source_errors, source_warnings = check_sources_tsv(ROOT, args.require_release_assets)
    errors.extend(source_errors)
    warnings.extend(source_warnings)

    selected_rows = [row for row in rows if row["status"] != "exclude"]
    final_rows = [row for row in selected_rows if row["status"] == "final"]
    provisional_rows = [
        row
        for row in selected_rows
        if row["status"] in {"provisional", "provisional-final"}
    ]

    print(f"Manifest rows: {len(rows)}")
    print(f"Selected rows: {len(selected_rows)}")
    print(f"Final rows: {len(final_rows)}")
    print(f"Provisional rows: {len(provisional_rows)}")
    print(f"Warnings: {len(warnings)}")
    print(f"Errors: {len(errors)}")

    for message in warnings:
        print(f"WARNING: {message}")
    for message in errors:
        print(f"ERROR: {message}")

    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
