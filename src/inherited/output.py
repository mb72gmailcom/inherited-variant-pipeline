from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

from inherited.constants import DEFAULT_BLOCK_SIZE


def serialize_patient_ids(hits: dict[str, tuple[str, str, str, str]]) -> str:
    """Serialize affected patient IDs only."""
    return ";".join(sorted(hits))


def serialize_trio_calls(
    hits: dict[str, tuple[str, str, str, str]],
) -> str:
    """Serialize {person_id: (mGT, fGT, cGT, cGQ)} to a single column string."""
    parts = []
    for person_id in sorted(hits):
        mother_gt, father_gt, child_gt, child_gq = hits[person_id]
        parts.append(f"{person_id}={mother_gt}|{father_gt}|{child_gt}|{child_gq}")
    return ";".join(parts)


def parse_trio_calls(payload: str) -> dict[str, tuple[str, str, str, str]]:
    hits: dict[str, tuple[str, str, str, str]] = {}
    if not payload:
        return hits
    for part in payload.split(";"):
        person_id, genotypes = part.split("=", 1)
        mother_gt, father_gt, child_gt, child_gq = genotypes.split("|")
        hits[person_id] = (mother_gt, father_gt, child_gt, child_gq)
    return hits


def serialize_payload(
    hits: dict[str, tuple[str, str, str, str]],
    *,
    short_format: bool,
) -> str:
    if short_format:
        return serialize_patient_ids(hits)
    return serialize_trio_calls(hits)


def parse_patient_ids(payload: str) -> list[str]:
    if not payload:
        return []
    return payload.split(";")


class BlockWriter:
    """Buffer tab-separated result lines and flush to disk in blocks."""

    HEADER_SHORT = "#CHROM\tPOS\tID\tREF\tALT\tPATIENTS\n"
    HEADER_FULL = "#CHROM\tPOS\tID\tREF\tALT\tTRIO_CALLS\n"

    def __init__(
        self,
        path: Path,
        block_size: int = DEFAULT_BLOCK_SIZE,
        *,
        short_format: bool = True,
    ) -> None:
        self.path = path
        self.block_size = block_size
        self.block: list[str] = []
        self.lines_written = 0
        self._handle = path.open("w", encoding="utf-8")
        header = self.HEADER_SHORT if short_format else self.HEADER_FULL
        self._handle.write(header)

    def append(self, chrom: str, pos: str, ref: str, alt: str, payload: str) -> None:
        self.block.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t{payload}")
        if len(self.block) >= self.block_size:
            self.flush()

    def flush(self) -> None:
        if not self.block:
            return
        self._handle.write("\n".join(self.block) + "\n")
        self.lines_written += len(self.block)
        self.block.clear()

    def close(self) -> None:
        self.flush()
        self._handle.close()


@dataclass
class ResultWriter:
    output_dir: Path
    block_size: int = DEFAULT_BLOCK_SIZE
    short_format: bool = True
    inherited_per_variant: dict[str, int] = field(default_factory=dict)
    inherited_per_person: dict[str, int] = field(default_factory=dict)
    mendelian_bad_per_gt: dict[str, int] = field(default_factory=dict)
    inherited_variants: int = 0
    mendelian_bad_variants: int = 0
    inherited_entries: int = 0
    mendelian_bad_entries: int = 0

    def __post_init__(self) -> None:
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.inherited = BlockWriter(
            self.output_dir / "inherited.tsv",
            self.block_size,
            short_format=self.short_format,
        )
        self.mendelian_bad = BlockWriter(
            self.output_dir / "mendelian_bad.tsv",
            self.block_size,
            short_format=self.short_format,
        )

    def write_inherited(
        self,
        chrom: str,
        pos: str,
        ref: str,
        alt: str,
        variant_key: str,
        hits: dict[str, tuple[str, str, str, str]],
    ) -> None:
        if not hits:
            return
        self.inherited.append(
            chrom, pos, ref, alt, serialize_payload(hits, short_format=self.short_format)
        )
        self.inherited_variants += 1
        self.inherited_entries += len(hits)
        self.inherited_per_variant[variant_key] = len(hits)
        for person_id in hits:
            self.inherited_per_person[person_id] = (
                self.inherited_per_person.get(person_id, 0) + 1
            )

    def write_mendelian_bad(
        self,
        chrom: str,
        pos: str,
        ref: str,
        alt: str,
        hits: dict[str, tuple[str, str, str, str]],
    ) -> None:
        if not hits:
            return
        self.mendelian_bad.append(
            chrom, pos, ref, alt, serialize_payload(hits, short_format=self.short_format)
        )
        self.mendelian_bad_variants += 1
        self.mendelian_bad_entries += len(hits)
        for _person_id, (mother_gt, father_gt, child_gt, _child_gq) in hits.items():
            gt_key = f"{mother_gt}:{father_gt}:{child_gt}"
            self.mendelian_bad_per_gt[gt_key] = self.mendelian_bad_per_gt.get(gt_key, 0) + 1

    def close(self) -> None:
        self.inherited.close()
        self.mendelian_bad.close()

    def save_summary_files(
        self,
        *,
        variants_seen: int,
        alleles_tested: int,
    ) -> None:
        _write_json(self.output_dir / "inherited_per_variant.json", self.inherited_per_variant)
        _write_json(self.output_dir / "inherited_per_person.json", self.inherited_per_person)
        _write_json(self.output_dir / "mendelian_bad_per_gt.json", self.mendelian_bad_per_gt)
        _write_json(
            self.output_dir / "stats.json",
            {
                "variants_seen": variants_seen,
                "alleles_tested": alleles_tested,
                "inherited_entries": self.inherited_entries,
                "inherited_variants": self.inherited_variants,
                "mendelian_bad_entries": self.mendelian_bad_entries,
                "mendelian_bad_variants": self.mendelian_bad_variants,
            },
        )


def _write_json(path: Path, payload: object) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)


def read_result_tsv(
    path: Path,
    *,
    short_format: bool = True,
) -> list[tuple[str, str, str, str, object]]:
    """Parse a result TSV into (chrom, pos, ref, alt, payload) records."""
    records = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            chrom, pos, _id, ref, alt, payload = line.rstrip("\n").split("\t", 5)
            if short_format:
                records.append((chrom, pos, ref, alt, parse_patient_ids(payload)))
            else:
                records.append((chrom, pos, ref, alt, parse_trio_calls(payload)))
    return records
