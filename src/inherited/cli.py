from __future__ import annotations

import argparse
import sys
from pathlib import Path

from inherited.analyze import analyze_vcf, save_run_params
from inherited.constants import DEFAULT_AF_THRESHOLD, DEFAULT_BLOCK_SIZE, DEFAULT_MEMORY_BLOCK


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="inherited",
        description="Classify inherited and mendelian-inconsistent rare variants in family trios",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    analyze = subparsers.add_parser(
        "analyze",
        help="Analyze a VCF using gnomAD AF JSON and family relations",
    )
    analyze.add_argument("--vcf", required=True, type=Path, help="Input VCF (.vcf or .vcf.gz)")
    analyze.add_argument(
        "--af-json",
        required=True,
        type=Path,
        help="JSON file with precomputed gnomAD allele frequencies",
    )
    analyze.add_argument(
        "--family-file",
        required=True,
        type=Path,
        help="Tab-separated family relations file",
    )
    analyze.add_argument(
        "-o",
        "--output-dir",
        required=True,
        type=Path,
        help="Directory for inherited.tsv and mendelian_bad.tsv",
    )
    analyze.add_argument(
        "--multiallelic",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Process multiallelic sites per alternate allele (default: True)",
    )
    analyze.add_argument(
        "--af-threshold",
        type=float,
        default=DEFAULT_AF_THRESHOLD,
        help=f"Maximum gnomAD AF to include (default: {DEFAULT_AF_THRESHOLD})",
    )
    analyze.add_argument(
        "--debug",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Print memory usage periodically during analysis (default: False)",
    )
    analyze.add_argument(
        "--memory-block",
        type=int,
        default=DEFAULT_MEMORY_BLOCK,
        help=f"Variant interval for debug memory logging (default: {DEFAULT_MEMORY_BLOCK})",
    )
    analyze.add_argument(
        "--block-size",
        type=int,
        default=DEFAULT_BLOCK_SIZE,
        help=f"Lines per block when streaming TSV output (default: {DEFAULT_BLOCK_SIZE})",
    )

    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)

    if args.command == "analyze":
        if not args.vcf.is_file():
            print(f"error: VCF not found: {args.vcf}", file=sys.stderr)
            raise SystemExit(1)
        if not args.af_json.is_file():
            print(f"error: AF JSON not found: {args.af_json}", file=sys.stderr)
            raise SystemExit(1)
        if not args.family_file.is_file():
            print(f"error: family file not found: {args.family_file}", file=sys.stderr)
            raise SystemExit(1)

        stats = analyze_vcf(
            vcf_path=args.vcf,
            af_json_path=args.af_json,
            family_file=args.family_file,
            output_dir=args.output_dir,
            multiallelic=args.multiallelic,
            af_threshold=args.af_threshold,
            debug=args.debug,
            memory_block=args.memory_block,
            block_size=args.block_size,
        )
        params_path = save_run_params(
            args.output_dir,
            vcf_path=args.vcf,
            af_json_path=args.af_json,
            family_file=args.family_file,
            multiallelic=args.multiallelic,
            af_threshold=args.af_threshold,
            debug=args.debug,
            memory_block=args.memory_block,
            block_size=args.block_size,
        )
        print(
            f"Wrote {stats.inherited_entries} inherited entries "
            f"({stats.inherited_variants} variants in inherited.tsv) and "
            f"{stats.mendelian_bad_entries} mendelian_bad entries "
            f"({stats.mendelian_bad_variants} variants in mendelian_bad.tsv) "
            f"to {args.output_dir}"
        )
        print(f"Wrote parameters to {params_path}")


if __name__ == "__main__":
    main()
