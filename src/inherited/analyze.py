from __future__ import annotations

import gzip
import json
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from inherited import __version__
from inherited.af import is_rare, load_af_json
from inherited.constants import (
    DEFAULT_AB,
    DEFAULT_AF_THRESHOLD,
    DEFAULT_DP,
    DEFAULT_GQ,
    DEFAULT_HAPLO_AB,
    DEFAULT_HAPLO_DP,
)
from inherited.families import build_trio_indices, load_family_relations
from inherited.genotype import get_good_site


def get_nfields(line: str, n: int) -> list[str]:
    return line.rstrip().split("\t")[:n]


@dataclass
class AnalysisStats:
    variants_seen: int = 0
    alleles_tested: int = 0
    inherited_entries: int = 0
    mendelian_bad_entries: int = 0
    inherited_variants: int = 0
    mendelian_bad_variants: int = 0


def summarize_inherited(
    dinh: dict[str, dict[str, tuple[str, str, str, str]]],
) -> tuple[dict[str, int], dict[str, int]]:
    """Return per-variant and per-person counts for inherited calls."""
    per_variant = {variant_key: len(people) for variant_key, people in dinh.items()}
    per_person: dict[str, int] = defaultdict(int)
    for people in dinh.values():
        for person_id in people:
            per_person[person_id] += 1
    return per_variant, dict(per_person)


def summarize_mendelian_bad(
    dm_bad: dict[str, dict[str, tuple[str, str, str, str]]],
) -> dict[str, int]:
    """Return counts keyed by mother_gt:father_gt:child_gt."""
    per_gt_pattern: dict[str, int] = defaultdict(int)
    for people in dm_bad.values():
        for _person_id, (mother_gt, father_gt, child_gt, _child_gq) in people.items():
            gt_key = f"{mother_gt}:{father_gt}:{child_gt}"
            per_gt_pattern[gt_key] += 1
    return dict(per_gt_pattern)


def analyze_vcf(
    vcf_path: Path,
    af_json_path: Path,
    family_file: Path,
    *,
    multiallelic: bool = True,
    af_threshold: float = DEFAULT_AF_THRESHOLD,
) -> tuple[dict[str, dict[str, tuple[str, str, str, str]]], dict[str, dict[str, tuple[str, str, str, str]]], AnalysisStats]:
    """Scan a VCF and classify rare variant trios into inherited and mendelian-bad buckets."""
    af_table = load_af_json(af_json_path)
    relations = load_family_relations(family_file)

    dinh: dict[str, dict[str, tuple[str, str, str, str]]] = defaultdict(dict)
    dm_bad: dict[str, dict[str, tuple[str, str, str, str]]] = defaultdict(dict)
    stats = AnalysisStats()

    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    with opener(vcf_path, "rt", encoding="utf-8") as handle:
        trios_ind: list[tuple[int, int, int]] = []
        sample_header: list[str] = []

        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                sample_header = line.strip().split("\t")[9:]
                _, trios_ind = build_trio_indices(sample_header, relations.trio_cl)
                continue

            if multiallelic:
                _process_multiallelic_line(
                    line,
                    af_table,
                    af_threshold,
                    trios_ind,
                    sample_header,
                    dinh,
                    dm_bad,
                    stats,
                )
            else:
                _process_biallelic_line(
                    line,
                    af_table,
                    af_threshold,
                    trios_ind,
                    sample_header,
                    dinh,
                    dm_bad,
                    stats,
                )

    return dict(dinh), dict(dm_bad), stats


def _process_multiallelic_line(
    line: str,
    af_table: dict[str, float],
    af_threshold: float,
    trios_ind: list[tuple[int, int, int]],
    sample_header: list[str],
    dinh: dict[str, dict[str, tuple[str, str, str, str]]],
    dm_bad: dict[str, dict[str, tuple[str, str, str, str]]],
    stats: AnalysisStats,
) -> None:
    chrom, pos, keys, ref, alts = get_nfields(line, 5)
    if len(ref) > 1:
        return

    skeys = keys.split(";")
    salts = alts.split(",")
    sample_fields = line.rstrip().split("\t")[9:]
    stats.variants_seen += 1

    for alt_index, key in enumerate(skeys, start=1):
        if alt_index > len(salts):
            break
        if not is_rare(af_table, key, af_threshold):
            continue

        alt = salts[alt_index - 1]
        if len(alt) > 1:
            continue

        stats.alleles_tested += 1
        _process_trios_for_allele(
            key,
            alt_index,
            sample_fields,
            sample_header,
            trios_ind,
            dinh,
            dm_bad,
            stats,
            clean_ad=True,
        )


def _process_biallelic_line(
    line: str,
    af_table: dict[str, float],
    af_threshold: float,
    trios_ind: list[tuple[int, int, int]],
    sample_header: list[str],
    dinh: dict[str, dict[str, tuple[str, str, str, str]]],
    dm_bad: dict[str, dict[str, tuple[str, str, str, str]]],
    stats: AnalysisStats,
) -> None:
    chrom, pos, key, ref, alt = get_nfields(line, 5)
    if len(ref) > 1 and len(alt) > 1:
        return
    if not is_rare(af_table, key, af_threshold):
        return

    stats.variants_seen += 1
    stats.alleles_tested += 1
    sample_fields = line.rstrip().split("\t")[9:]
    _process_trios_for_allele(
        key,
        1,
        sample_fields,
        sample_header,
        trios_ind,
        dinh,
        dm_bad,
        stats,
    )


def _process_trios_for_allele(
    key: str,
    alt_index: int,
    sample_fields: list[str],
    sample_header: list[str],
    trios_ind: list[tuple[int, int, int]],
    dinh: dict[str, dict[str, tuple[str, str, str, str]]],
    dm_bad: dict[str, dict[str, tuple[str, str, str, str]]],
    stats: AnalysisStats,
    *,
    clean_ad: bool = False,
) -> None:
    parents_cache: dict[int, list[object]] = {}

    for child_idx, mother_idx, father_idx in trios_ind:
        child_sample = sample_fields[child_idx]
        ac, child_gt, child_gq = get_good_site(child_sample, alt_index, clean_ad=clean_ad)
        if ac <= 0:
            continue

        if mother_idx in parents_cache:
            mac, mother_gt, mother_gq = parents_cache[mother_idx]
        else:
            mac, mother_gt, mother_gq = get_good_site(
                sample_fields[mother_idx], alt_index, clean_ad=clean_ad
            )
            parents_cache[mother_idx] = [mac, mother_gt, mother_gq]

        if father_idx in parents_cache:
            fac, father_gt, father_gq = parents_cache[father_idx]
        else:
            fac, father_gt, father_gq = get_good_site(
                sample_fields[father_idx], alt_index, clean_ad=clean_ad
            )
            parents_cache[father_idx] = [fac, father_gt, father_gq]

        pid = sample_header[child_idx]
        bucket: str | None = None
        record: tuple[str, str, str, str]

        if mac > 0 and fac > 0:
            if mac == 2 and fac == 2 and ac < 2:
                bucket = "mendelian_bad"
            else:
                bucket = "inherited"
            record = (mother_gt, father_gt, child_gt, child_gq)
        elif mac > 0:
            bucket = "mendelian_bad" if ac == 2 else "inherited"
            record = (mother_gt, "0/0", child_gt, child_gq)
        elif fac > 0:
            bucket = "mendelian_bad" if ac == 2 else "inherited"
            record = ("0/0", father_gt, child_gt, child_gq)
        else:
            continue

        if bucket == "inherited":
            dinh[key][pid] = record
            stats.inherited_entries += 1
        else:
            dm_bad[key][pid] = record
            stats.mendelian_bad_entries += 1


def save_results(
    output_dir: Path,
    dinh: dict[str, dict[str, tuple[str, str, str, str]]],
    dm_bad: dict[str, dict[str, tuple[str, str, str, str]]],
    stats: AnalysisStats,
) -> AnalysisStats:
    output_dir.mkdir(parents=True, exist_ok=True)

    inherited_per_variant, inherited_per_person = summarize_inherited(dinh)
    mendelian_bad_per_gt = summarize_mendelian_bad(dm_bad)

    stats.inherited_variants = len(dinh)
    stats.mendelian_bad_variants = len(dm_bad)

    _write_json(output_dir / "inherited.json", dinh)
    _write_json(output_dir / "mendelian_bad.json", dm_bad)
    _write_json(output_dir / "inherited_per_variant.json", inherited_per_variant)
    _write_json(output_dir / "inherited_per_person.json", inherited_per_person)
    _write_json(output_dir / "mendelian_bad_per_gt.json", mendelian_bad_per_gt)
    _write_json(
        output_dir / "stats.json",
        {
            "variants_seen": stats.variants_seen,
            "alleles_tested": stats.alleles_tested,
            "inherited_entries": stats.inherited_entries,
            "inherited_variants": stats.inherited_variants,
            "mendelian_bad_entries": stats.mendelian_bad_entries,
            "mendelian_bad_variants": stats.mendelian_bad_variants,
        },
    )
    return stats


def save_run_params(
    output_dir: Path,
    *,
    vcf_path: Path,
    af_json_path: Path,
    family_file: Path,
    multiallelic: bool,
    af_threshold: float,
) -> Path:
    """Write the parameters for this run into the chromosome output directory."""
    output_dir.mkdir(parents=True, exist_ok=True)
    params_path = output_dir / "params.json"

    payload: dict[str, Any] = {
        "package_version": __version__,
        "run_at": datetime.now(timezone.utc).isoformat(),
        "vcf": str(vcf_path.resolve()),
        "af_json": str(af_json_path.resolve()),
        "family_file": str(family_file.resolve()),
        "output_dir": str(output_dir.resolve()),
        "multiallelic": multiallelic,
        "af_threshold": af_threshold,
        "quality_filters": {
            "gq": DEFAULT_GQ,
            "dp": DEFAULT_DP,
            "ab": DEFAULT_AB,
            "haplo_dp": DEFAULT_HAPLO_DP,
            "haplo_ab": DEFAULT_HAPLO_AB,
        },
    }
    _write_json(params_path, payload)
    return params_path


def _write_json(path: Path, payload: object) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
