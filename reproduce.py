#!/usr/bin/env python3
"""Reproduce released SORT1 analyses from AlphaGenome and public raw inputs."""

from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

from reproduction.common import (
    Audit,
    REPO_ROOT,
    archive_previous_audit,
    credential_present,
    initialize_run,
    load_api_key_file,
    load_env_file,
    make_alphagenome_importable,
)
from reproduction.figure1 import (
    build_fig1c_variant_table,
    fetch_fig1c_inputs,
    run_fig1b,
    run_fig1c_middle,
    run_fig1e,
    validate_fig1c_reference,
)
from reproduction.figure1_public import (
    build_full_figure1c,
    prepare_figure1_public_inputs,
    run_observed_contact,
)
from reproduction.figure2 import (
    fetch_hg38,
    fetch_kircher,
    prepare_figure2_inputs,
    run_fig2b,
    run_fig2c,
    run_fig2e,
    run_fig2f,
)
from reproduction.figure3 import run_fig3a, run_fig3b, run_fig3c, run_fig3e, run_fig3f, run_fig3g
from reproduction.figure4 import run_fig4b, run_fig4c
from reproduction.figure4ef import run_fig4ef
from reproduction.figure4g import run_fig4g
from reproduction.figure4h import run_fig4h
from reproduction.figureS1 import run_figs1a, run_figs1b, run_figs1c, run_figs1d
from reproduction.figureS2 import run_figs2a, run_figs2b, run_figs2c
from reproduction.figureS3 import run_figs3
from reproduction.figureS4 import run_figs4a, run_figs4b, run_figs4c, run_figs4d
from reproduction.figureS5 import run_figs5_substitutions, run_figs5d
from reproduction.report import compare_run, write_report


FIGURE1_PANELS = ["1B", "1C", "1D", "1E", "1F"]
FIGURE2_PANELS = ["2B", "2C", "2E", "2F"]
FIGURE3_PANELS = ["3A", "3B", "3C", "3E", "3F", "3G"]
FIGURE4_PANELS = ["4B", "4C", "4E", "4F", "4G", "4H"]
SUPP_S1_PANELS = ["S1A", "S1B", "S1C", "S1D"]
SUPP_S2_PANELS = ["S2A", "S2B", "S2C"]
SUPP_S3_PANELS = ["S3A", "S3B", "S3C", "S3D", "S3E", "S3F", "S3G"]
SUPP_S4_PANELS = ["S4A", "S4B", "S4C", "S4D"]
SUPP_S5_PANELS = ["S5A", "S5B", "S5C", "S5D"]
SUPPORTED_PANELS = (
    FIGURE1_PANELS + ["1C-middle"] + FIGURE2_PANELS + FIGURE3_PANELS + FIGURE4_PANELS
    + SUPP_S1_PANELS + SUPP_S2_PANELS + SUPP_S3_PANELS + SUPP_S4_PANELS + SUPP_S5_PANELS
)
DEFAULT_PANELS = FIGURE1_PANELS
DEFAULT_KEY_FILE = REPO_ROOT / "ALPHAGENOME_API_KEY.txt"


def parse_panels(value: str) -> list[str]:
    panels = [item.strip() for item in value.split(",") if item.strip()]
    unknown = sorted(set(panels) - set(SUPPORTED_PANELS))
    if unknown:
        raise argparse.ArgumentTypeError(f"Unsupported panel(s): {', '.join(unknown)}")
    return panels


def default_run_dir(panels: list[str] | None = None) -> Path:
    figure = (
        "figure4" if panels and set(panels).issubset(FIGURE4_PANELS)
        else "figure3" if panels and set(panels).issubset(FIGURE3_PANELS)
        else "figure2" if panels and set(panels).issubset(FIGURE2_PANELS)
        else "supp_s1" if panels and set(panels).issubset(SUPP_S1_PANELS)
        else "supp_s2" if panels and set(panels).issubset(SUPP_S2_PANELS)
        else "supp_s3" if panels and set(panels).issubset(SUPP_S3_PANELS)
        else "supp_s4" if panels and set(panels).issubset(SUPP_S4_PANELS)
        else "supp_s5" if panels and set(panels).issubset(SUPP_S5_PANELS)
        else "reproduction"
    )
    stamp = datetime.now(timezone.utc).strftime(f"{figure}_%Y%m%dT%H%M%SZ")
    return REPO_ROOT / "reproduction_runs" / stamp


def doctor(env_file: Path, api_key_file: Path | None) -> int:
    loaded = load_env_file(env_file)
    key_file_loaded = load_api_key_file(api_key_file)
    source = make_alphagenome_importable()
    checks: dict[str, dict[str, object]] = {}
    for module in ("numpy", "pandas", "matplotlib", "pysam", "pyliftover", "alphagenome"):
        try:
            imported = __import__(module)
            checks[module] = {"ok": True, "version": getattr(imported, "__version__", "installed")}
        except Exception as exc:
            checks[module] = {"ok": False, "error": f"{type(exc).__name__}: {exc}"}
    checks["credential"] = {
        "ok": credential_present(),
        "source": "environment or .env (presence only; not service-validated)"
        if credential_present()
        else "missing or placeholder",
    }
    result = {
        "ok": all(bool(item["ok"]) for item in checks.values()),
        "repository": str(REPO_ROOT),
        "env_file": str(env_file),
        "loaded_names": loaded,
        "api_key_file_loaded": key_file_loaded,
        "development_sdk_source": str(source) if source else None,
        "checks": checks,
    }
    print(json.dumps(result, indent=2))
    return 0 if result["ok"] else 1


def run(args: argparse.Namespace) -> int:
    load_env_file(args.env_file)
    load_api_key_file(args.api_key_file)
    initialize_run(args.run_dir, resume=args.resume)
    if args.resume:
        archive_previous_audit(args.run_dir)
    audit = Audit(args.run_dir.resolve(), args.panels)
    audit.save()
    try:
        figure1_inputs: dict[str, Path] = {}
        if "1C" in args.panels:
            with audit.step("Figure 1C: stage GTEx v7 and 1000 Genomes public inputs"):
                figure1_inputs = prepare_figure1_public_inputs(
                    args.run_dir, audit, gtex_file=args.gtex_file,
                    vcf_file=args.onekg_vcf, panel_file=args.onekg_panel,
                )
        figure2_inputs: dict[str, Path] = {}
        figure3_fasta: Path | None = None
        if set(args.panels) & set(FIGURE3_PANELS):
            with audit.step("Figure 3: download and checksum GRCh38"):
                figure3_fasta = fetch_hg38(args.run_dir, audit)
        if set(args.panels) & set(FIGURE2_PANELS):
            with audit.step("Figure 2: download and checksum original public inputs"):
                figure2_inputs = prepare_figure2_inputs(
                    args.run_dir, audit, args.panels, args.wang_xls
                )
        if "1B" in args.panels:
            with audit.step("1B: AlphaGenome locus tracks and rendering"):
                run_fig1b(args.run_dir, audit)
        if "1C-middle" in args.panels or "1C" in args.panels:
            run_fig1c_middle(
                args.run_dir,
                audit,
                batch_size=args.batch_size,
                max_workers=args.max_workers,
                max_variants=args.max_variants,
            )
        if "1C" in args.panels:
            with audit.step("1C: reconstruct GTEx and EUR-tagging panels"):
                build_full_figure1c(args.run_dir, figure1_inputs)
        if "1E" in args.panels:
            with audit.step("1E: AlphaGenome FOLD_0 contact map and rendering"):
                run_fig1e(args.run_dir, audit)
        if "1D" in args.panels or "1F" in args.panels:
            if "1F" in args.panels and not (args.run_dir / "predictions/Figure1E_fold0_contact/ref_alt.npz").exists():
                with audit.step("1F prerequisite: AlphaGenome FOLD_0 contact map"):
                    run_fig1e(args.run_dir, audit)
            with audit.step("1D/1F: extract public 4DN HepG2 observed Hi-C"):
                run_observed_contact(args.run_dir, audit, args.hic_file)
        if "2B" in args.panels:
            with audit.step("2B: reconstruct Wang repair products and score liver RNA"):
                run_fig2b(
                    args.run_dir,
                    audit,
                    figure2_inputs["fasta"],
                    figure2_inputs["wang"],
                    batch_size=args.batch_size,
                    max_workers=args.max_workers,
                )
        if "2C" in args.panels:
            with audit.step("2C: construct and score the 11-by-11 deletion grid"):
                run_fig2c(
                    args.run_dir,
                    audit,
                    figure2_inputs["fasta"],
                    batch_size=args.batch_size,
                    max_workers=args.max_workers,
                )
        if "2E" in args.panels:
            with audit.step("2E: score the 1,798-substitution SORT1 Kircher benchmark"):
                run_fig2e(
                    args.run_dir,
                    audit,
                    figure2_inputs["kircher"],
                    chunk_size=args.ism_chunk_size,
                    max_workers=args.ism_max_workers,
                )
        if "2F" in args.panels:
            with audit.step("2F: score the six-element Kircher accessibility benchmark"):
                run_fig2f(args.run_dir, audit, figure2_inputs["kircher"])
        if "3B" in args.panels:
            with audit.step("3B: native-locus 501-bp three-gene ISM"):
                run_fig3b(args.run_dir, audit, figure3_fasta, batch_size=args.batch_size, max_workers=args.max_workers)
        if "3C" in args.panels:
            if not (args.run_dir / "derived/Figure3B_native_501bp_ism/native_locus_501bp_all_gene_scores.tsv").exists():
                with audit.step("3C prerequisite: Figure 3B native-locus 501-bp three-gene ISM"):
                    run_fig3b(args.run_dir, audit, figure3_fasta, batch_size=args.batch_size, max_workers=args.max_workers)
            with audit.step("3C: scan JASPAR motifs and score PWM compatibility vs Figure 3B ISM"):
                run_fig3c(args.run_dir, audit, figure3_fasta)
        if "3G" in args.panels:
            with audit.step("3G: expanded component-necessity audit"):
                run_fig3g(args.run_dir, audit, figure3_fasta, batch_size=args.batch_size, max_workers=args.max_workers)
        if "3E" in args.panels:
            with audit.step("3E: directional single-arm motif-protected recovery"):
                run_fig3e(args.run_dir, audit, figure3_fasta, batch_size=args.batch_size, max_workers=args.max_workers)
        if "3A" in args.panels:
            with audit.step("3A: 100kb regional two-stage RNA(TSS) ISM scan"):
                run_fig3a(args.run_dir, audit, figure3_fasta)
        if "3F" in args.panels:
            with audit.step("3F: wide-main-panel 1bp boundary grid"):
                run_fig3f(args.run_dir, audit, figure3_fasta, batch_size=args.batch_size, max_workers=args.max_workers)
        if "4B" in args.panels:
            with audit.step("4B: bottom100 315bp donor eight-distance sweep"):
                run_fig4b(args.run_dir, audit, batch_size=args.batch_size, max_workers=args.max_workers)
        if "4C" in args.panels:
            with audit.step("4C: HPA bottom/middle/top-500 cohort fold-change at 30bp"):
                run_fig4c(args.run_dir, audit, batch_size=args.batch_size, max_workers=args.max_workers, hpa_file=args.hpa_file)
        if "4E" in args.panels or "4F" in args.panels:
            done_marker = args.run_dir / "derived/Figure4EF_distal_contact_transfer/Figure4F_contact_dose_response.tsv"
            if not done_marker.exists():
                with audit.step("4E/4F: chr1 Hi-C-guided distal-contact 315bp T/G transfer benchmark"):
                    run_fig4ef(args.run_dir, audit, batch_size=args.batch_size, max_workers=args.max_workers)
        if "4G" in args.panels:
            run_fig4g(args.run_dir, audit)
        if "4H" in args.panels:
            with audit.step("4H: exhaustive +/-50kb x 3-alt x 3-tissue regional ISM synergy scan"):
                run_fig4h(args.run_dir, audit, batch_size=args.batch_size, max_workers=args.max_workers)
        if "S1A" in args.panels:
            run_figs1a(args.run_dir, audit)
        if "S1B" in args.panels:
            run_figs1b(args.run_dir, audit)
        if "S1C" in args.panels:
            run_figs1c(args.run_dir, audit)
        if "S1D" in args.panels:
            run_figs1d(args.run_dir, audit)
        if "S2A" in args.panels:
            run_figs2a(args.run_dir, audit)
        if "S2B" in args.panels:
            run_figs2b(args.run_dir, audit)
        if "S2C" in args.panels:
            run_figs2c(args.run_dir, audit)
        if set(args.panels) & set(SUPP_S3_PANELS):
            run_figs3(args.run_dir, audit, currin_variants_file=args.currin_variants, currin_peakset_file=args.currin_peakset)
        if "S4A" in args.panels:
            run_figs4a(args.run_dir, audit, wang_xls=args.wang_xls)
        if "S4B" in args.panels:
            run_figs4b(args.run_dir, audit, wang_xls=args.wang_xls)
        if "S4C" in args.panels:
            run_figs4c(args.run_dir, audit)
        if "S4D" in args.panels:
            run_figs4d(args.run_dir, audit)
        if set(args.panels) & {"S5A", "S5B", "S5C"}:
            kircher_path = figure2_inputs.get("kircher") or fetch_kircher(args.run_dir, audit)
            run_figs5_substitutions(args.run_dir, audit, kircher_path, chunk_size=args.ism_chunk_size, max_workers=args.ism_max_workers)
        if "S5D" in args.panels:
            kircher_path = figure2_inputs.get("kircher") or fetch_kircher(args.run_dir, audit)
            fasta_path = figure2_inputs.get("fasta") or fetch_hg38(args.run_dir, audit)
            run_figs5d(args.run_dir, audit, kircher_path, fasta_path)
        audit.finish()
    except BaseException:
        write_report(args.run_dir)
        raise
    report = write_report(args.run_dir)
    print(f"Run complete: {args.run_dir}")
    print(f"Audit report: {report}")
    return 0


def prepare(args: argparse.Namespace) -> int:
    """Download and reconstruct public inputs without requiring an API key."""
    initialize_run(args.run_dir, resume=args.resume)
    if args.resume:
        archive_previous_audit(args.run_dir)
    audit = Audit(args.run_dir.resolve(), args.panels)
    audit.save()
    try:
        if "1C-middle" in args.panels or "1C" in args.panels:
            with audit.step("1C-middle: download GLGC, LiftOver chain, and hg38"):
                inputs = fetch_fig1c_inputs(args.run_dir, audit)
            with audit.step("1C-middle: reconstruct the 111-variant set") as record:
                variants = build_fig1c_variant_table(
                    inputs["glgc"],
                    inputs["chain"],
                    args.run_dir / "derived" / "Figure1C_reconstructed_variants.tsv",
                )
                record["variants"] = int(len(variants))
                record["grch38_reference_allele_matches"] = validate_fig1c_reference(
                    variants, inputs["fasta"]
                )
        if "1C" in args.panels:
            with audit.step("Figure 1C: stage GTEx v7 and 1000 Genomes public inputs"):
                prepare_figure1_public_inputs(
                    args.run_dir, audit, gtex_file=args.gtex_file,
                    vcf_file=args.onekg_vcf, panel_file=args.onekg_panel,
                )
        if set(args.panels) & set(FIGURE2_PANELS):
            with audit.step("Figure 2: download and checksum original public inputs"):
                prepare_figure2_inputs(args.run_dir, audit, args.panels, args.wang_xls)
        if set(args.panels) & set(FIGURE3_PANELS):
            with audit.step("Figure 3: download and checksum GRCh38"):
                fetch_hg38(args.run_dir, audit)
        audit.finish()
    except BaseException:
        write_report(args.run_dir)
        raise
    report = write_report(args.run_dir)
    print(f"Public inputs prepared: {args.run_dir}")
    print(f"Audit report: {report}")
    return 0


def compare(args: argparse.Namespace) -> int:
    run_json = args.run_dir / "audit" / "run.json"
    if not run_json.exists():
        raise FileNotFoundError(f"Not a reproduction run directory: {args.run_dir}")
    audit = json.loads(run_json.read_text())
    panels = args.panels or audit["panels"]
    result = compare_run(args.run_dir, panels, reference_4h_file=args.reference_4h_file)
    report = write_report(args.run_dir, result)
    print(f"Comparison: {'PASS' if result['pass'] else 'FAIL'}")
    print(f"Report: {report}")
    return 0 if result["pass"] else 2


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--env-file", type=Path, default=REPO_ROOT / ".env")
    parser.add_argument(
        "--api-key-file",
        type=Path,
        default=None,
        help="Optional one-line credential file; its value is never printed or copied.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("doctor", help="Check dependencies and credential presence")
    prepare_parser = subparsers.add_parser(
        "prepare", help="Download and validate public inputs without API access"
    )
    prepare_parser.add_argument("--panels", type=parse_panels, default=FIGURE1_PANELS)
    prepare_parser.add_argument("--run-dir", type=Path, default=None)
    prepare_parser.add_argument("--resume", action="store_true")
    prepare_parser.add_argument(
        "--wang-xls", type=Path, default=None,
        help="Manually downloaded original atv310103_ds.xls if the publisher blocks automation.",
    )
    run_parser = subparsers.add_parser("run", help="Create a fresh analysis run")
    run_parser.add_argument("--panels", type=parse_panels, default=DEFAULT_PANELS)
    run_parser.add_argument("--run-dir", type=Path, default=None)
    run_parser.add_argument("--resume", action="store_true")
    run_parser.add_argument("--batch-size", type=int, default=8)
    run_parser.add_argument("--max-workers", type=int, default=4)
    run_parser.add_argument("--ism-chunk-size", type=int, default=10)
    run_parser.add_argument("--ism-max-workers", type=int, default=1)
    run_parser.add_argument(
        "--wang-xls", type=Path, default=None,
        help="Manually downloaded original atv310103_ds.xls if the publisher blocks automation.",
    )
    run_parser.add_argument(
        "--hpa-file", type=Path, default=None,
        help="Manually downloaded HPA v24.1 rna_tissue_consensus.tsv if the "
             "portal blocks automation (checksum-verified against data/SOURCES.tsv).",
    )
    run_parser.add_argument(
        "--max-variants",
        type=int,
        default=None,
        help="Development smoke test only; a limited run cannot pass publication comparison.",
    )
    for panel_parser in (prepare_parser, run_parser):
        panel_parser.add_argument("--gtex-file", type=Path, default=None,
            help="Optional manually downloaded GTEx v7 Liver.allpairs.txt.gz (about 3.4 GB).")
        panel_parser.add_argument("--onekg-vcf", type=Path, default=None,
            help="Optional phased 1000 Genomes chr1/SORT1 interval VCF.gz; its .tbi is copied when present.")
        panel_parser.add_argument("--onekg-panel", type=Path, default=None,
            help="Optional Phase 3 1000 Genomes sample panel file.")
    run_parser.add_argument("--hic-file", type=Path, default=None,
        help="Optional local 4DNFICSTCJQZ.hic; default uses remote HTTP byte-range access.")
    run_parser.add_argument("--currin-variants", type=Path, default=None,
        help="Required for Figure S3: local Currin et al. 2025 SORT1-locus caQTL "
             "association table (Zenodo 15025748 / GEO GSE264684; checksum-verified).")
    run_parser.add_argument("--currin-peakset", type=Path, default=None,
        help="Required for Figure S3: local Currin et al. 2025 28-peak coordinated-set "
             "definition for the SORT1 locus (checksum-verified).")
    compare_parser = subparsers.add_parser("compare", help="Compare a completed run to frozen references")
    compare_parser.add_argument("--run-dir", type=Path, required=True)
    compare_parser.add_argument("--panels", type=parse_panels, default=None)
    compare_parser.add_argument(
        "--reference-4h-file", type=Path, default=None,
        help="Manually supplied copy of Figure4H_regional_tissue_scan.tsv (checksum-verified) "
             "since that reference is Zenodo-pending and not committed to this repository.",
    )
    args = parser.parse_args()
    if args.api_key_file is None and DEFAULT_KEY_FILE.exists():
        args.api_key_file = DEFAULT_KEY_FILE
    if args.command == "doctor":
        raise SystemExit(doctor(args.env_file, args.api_key_file))
    if args.command == "prepare":
        args.run_dir = (args.run_dir or default_run_dir(args.panels)).resolve()
        raise SystemExit(prepare(args))
    if args.command == "run":
        args.run_dir = (args.run_dir or default_run_dir(args.panels)).resolve()
        raise SystemExit(run(args))
    raise SystemExit(compare(args))


if __name__ == "__main__":
    main()
