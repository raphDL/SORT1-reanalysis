from __future__ import annotations

import os
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from reproduce import parse_panels
from reproduction.common import (
    archive_previous_audit,
    file_manifest,
    initialize_run,
    load_api_key_file,
    load_env_file,
)
from reproduction.figure2 import _minor_baseline
from reproduction.report import _numeric_summary


class ReproductionTests(unittest.TestCase):
    def test_panel_parser(self) -> None:
        self.assertEqual(parse_panels("1B,1C-middle,1E"), ["1B", "1C-middle", "1E"])
        self.assertEqual(parse_panels("2B,2C,2E,2F"), ["2B", "2C", "2E", "2F"])

    def test_run_directory_must_be_clean(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "run"
            initialize_run(root, resume=False)
            with self.assertRaises(FileExistsError):
                initialize_run(root, resume=False)
            initialize_run(root, resume=True)

    def test_env_loader_does_not_override_shell(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / ".env"
            path.write_text("SORT1_TEST_SECRET=from_file\n")
            os.environ["SORT1_TEST_SECRET"] = "from_shell"
            try:
                self.assertEqual(load_env_file(path), [])
                self.assertEqual(os.environ["SORT1_TEST_SECRET"], "from_shell")
            finally:
                os.environ.pop("SORT1_TEST_SECRET", None)

    def test_assignment_style_key_file(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "key.txt"
            path.write_text('ALPHAGENOME_API_KEY="test-value"\n')
            os.environ.pop("ALPHAGENOME_API_KEY", None)
            try:
                self.assertTrue(load_api_key_file(path))
                self.assertEqual(os.environ["ALPHAGENOME_API_KEY"], "test-value")
            finally:
                os.environ.pop("ALPHAGENOME_API_KEY", None)

    def test_numeric_comparison_tolerance(self) -> None:
        result = _numeric_summary(np.array([1.0, 2.0]), np.array([1.0, 2.0 + 1e-7]))
        self.assertTrue(result["pass"])

    def test_numeric_comparison_records_panel_specific_threshold(self) -> None:
        result = _numeric_summary(
            np.array([1.0, 2.1, 3.0]), np.array([1.0, 2.0, 3.0]),
            rtol=0.0, atol=0.11, min_pearson=0.99,
        )
        self.assertTrue(result["pass"])
        self.assertEqual(result["atol"], 0.11)
        self.assertEqual(result["minimum_pearson_r"], 0.99)

    def test_resumed_audit_is_archived(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "audit").mkdir()
            (root / "audit" / "run.json").write_text('{"started_utc":"2026-08-02T09:00:00+00:00"}\n')
            (root / "audit" / "REPRODUCIBILITY_REPORT.md").write_text("old report\n")
            (root / "audit" / "comparison.json").write_text('{"pass":true}\n')
            archive_previous_audit(root)
            archived = list((root / "audit" / "attempts").glob("*/run.json"))
            self.assertEqual(len(archived), 1)
            self.assertFalse((root / "audit" / "comparison.json").exists())
            self.assertEqual(len(list((root / "audit" / "attempts").glob("*/comparison.json"))), 1)

    def test_fig2c_minor_baseline_survives_duplicate_zero_deletion_row(self) -> None:
        # Regression test for the intermittent "Merge keys are not unique in
        # right dataset; not a many-to-one merge" MergeError: the grid cell
        # (upstream=0, downstream=0) deletes zero bases and hashes identically
        # to the standalone "minor" design, so two design rows legitimately
        # share `minor_hash` and every `keys` value.
        minor_hash = "minor_hash_value"
        keys = ["gene_symbol", "target_index"]
        expanded = pd.DataFrame(
            [
                {"sequence_sha256": minor_hash, "gene_symbol": "SORT1", "target_index": 4,
                 "design_id": "minor", "rna_mean_tss_pm2kb": 0.5},
                {"sequence_sha256": minor_hash, "gene_symbol": "SORT1", "target_index": 4,
                 "design_id": "del_xy_u00_d00", "rna_mean_tss_pm2kb": 0.5},
                {"sequence_sha256": "other_hash", "gene_symbol": "SORT1", "target_index": 4,
                 "design_id": "del_xy_u01_d00", "rna_mean_tss_pm2kb": 0.6},
            ]
        )
        baseline = _minor_baseline(expanded, keys, minor_hash)
        self.assertEqual(len(baseline), 1)
        self.assertAlmostEqual(baseline.rna_mean_tss_pm2kb_minor.iloc[0], 0.5)
        # The result must merge cleanly as many_to_one against the full frame.
        expanded.merge(baseline, on=keys, how="left", validate="many_to_one")

    def test_fig2c_minor_baseline_survives_nan_optional_metadata_column(self) -> None:
        minor_hash = "minor_hash_value"
        keys = ["gene_symbol", "biosample_name"]
        expanded = pd.DataFrame(
            [
                {"sequence_sha256": minor_hash, "gene_symbol": "SORT1", "biosample_name": np.nan,
                 "rna_mean_tss_pm2kb": 0.5},
                {"sequence_sha256": minor_hash, "gene_symbol": "SORT1", "biosample_name": np.nan,
                 "rna_mean_tss_pm2kb": 0.5},
            ]
        )
        baseline = _minor_baseline(expanded, keys, minor_hash)
        self.assertEqual(len(baseline), 1)

    def test_fig2c_minor_baseline_raises_on_genuine_disagreement(self) -> None:
        minor_hash = "minor_hash_value"
        keys = ["gene_symbol"]
        expanded = pd.DataFrame(
            [
                {"sequence_sha256": minor_hash, "gene_symbol": "SORT1", "rna_mean_tss_pm2kb": 0.5},
                {"sequence_sha256": minor_hash, "gene_symbol": "SORT1", "rna_mean_tss_pm2kb": 0.9},
            ]
        )
        with self.assertRaises(ValueError):
            _minor_baseline(expanded, keys, minor_hash)

    def test_manifest_excludes_mutable_audit_files(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "audit").mkdir()
            (root / "audit" / "run.json").write_text("{}")
            (root / "audit" / "REPRODUCIBILITY_REPORT.md").write_text("report")
            (root / "result.tsv").write_text("value\n")
            self.assertEqual([item["path"] for item in file_manifest(root)], ["result.tsv"])


if __name__ == "__main__":
    unittest.main()
