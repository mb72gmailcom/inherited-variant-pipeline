from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

from inherited.checkpoint import (
    STATS_CUMULATIVE_FILENAME,
    Checkpoint,
    CumulativeStats,
    save_checkpoint,
    write_json_atomic,
)
from inherited.constants import DEFAULT_BLOCK_SIZE, DEFAULT_SEGMENT_SIZE


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
    segment_size: int = DEFAULT_SEGMENT_SIZE
    short_format: bool = True
    cumulative: CumulativeStats = field(default_factory=CumulativeStats)
    segment_index: int = 0
    last_chrom: str = ""
    last_pos: int = 0
    inherited_segment_lines: int = 0
    mendelian_bad_segment_lines: int = 0
    segment_inherited_per_variant: dict[str, int] = field(default_factory=dict)
    _inherited: BlockWriter | None = field(default=None, repr=False)
    _mendelian_bad: BlockWriter | None = field(default=None, repr=False)

    def __post_init__(self) -> None:
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self._open_segment_writers()

    @classmethod
    def from_checkpoint(
        cls,
        output_dir: Path,
        checkpoint: Checkpoint,
        *,
        block_size: int = DEFAULT_BLOCK_SIZE,
        segment_size: int = DEFAULT_SEGMENT_SIZE,
        short_format: bool = True,
    ) -> ResultWriter:
        writer = cls(
            output_dir=output_dir,
            block_size=block_size,
            segment_size=segment_size,
            short_format=short_format,
            cumulative=checkpoint.cumulative,
            segment_index=checkpoint.segment_index + 1,
            last_chrom=checkpoint.chrom,
            last_pos=checkpoint.last_pos,
        )
        return writer

    def _inherited_path(self, segment_index: int) -> Path:
        if self.segment_size <= 0:
            return self.output_dir / "inherited.tsv"
        return self.output_dir / f"inherited_{segment_index:05d}.tsv"

    def _mendelian_bad_path(self, segment_index: int) -> Path:
        if self.segment_size <= 0:
            return self.output_dir / "mendelian_bad.tsv"
        return self.output_dir / f"mendelian_bad_{segment_index:05d}.tsv"

    def _open_segment_writers(self) -> None:
        if self._inherited is not None:
            self._inherited.close()
        if self._mendelian_bad is not None:
            self._mendelian_bad.close()
        self._inherited = BlockWriter(
            self._inherited_path(self.segment_index),
            self.block_size,
            short_format=self.short_format,
        )
        self._mendelian_bad = BlockWriter(
            self._mendelian_bad_path(self.segment_index),
            self.block_size,
            short_format=self.short_format,
        )
        self.inherited_segment_lines = 0
        self.mendelian_bad_segment_lines = 0

    def write_inherited(
        self,
        chrom: str,
        pos: str,
        ref: str,
        alt: str,
        variant_key: str,
        hits: dict[str, tuple[str, str, str, str]],
    ) -> None:
        if not hits or self._inherited is None:
            return
        self._inherited.append(
            chrom, pos, ref, alt, serialize_payload(hits, short_format=self.short_format)
        )
        self.inherited_segment_lines += 1
        self.cumulative.inherited_variants += 1
        self.cumulative.inherited_entries += len(hits)
        self.segment_inherited_per_variant[variant_key] = len(hits)
        for person_id in hits:
            self.cumulative.inherited_per_person[person_id] = (
                self.cumulative.inherited_per_person.get(person_id, 0) + 1
            )
        self._update_position(chrom, pos)
        self._maybe_rotate_segment()

    def write_mendelian_bad(
        self,
        chrom: str,
        pos: str,
        ref: str,
        alt: str,
        hits: dict[str, tuple[str, str, str, str]],
    ) -> None:
        if not hits or self._mendelian_bad is None:
            return
        self._mendelian_bad.append(
            chrom, pos, ref, alt, serialize_payload(hits, short_format=self.short_format)
        )
        self.mendelian_bad_segment_lines += 1
        self.cumulative.mendelian_bad_variants += 1
        self.cumulative.mendelian_bad_entries += len(hits)
        for _person_id, (mother_gt, father_gt, child_gt, _child_gq) in hits.items():
            gt_key = f"{mother_gt}:{father_gt}:{child_gt}"
            self.cumulative.mendelian_bad_per_gt[gt_key] = (
                self.cumulative.mendelian_bad_per_gt.get(gt_key, 0) + 1
            )
        self._update_position(chrom, pos)
        self._maybe_rotate_segment()

    def _update_position(self, chrom: str, pos: str) -> None:
        self.last_chrom = chrom
        self.last_pos = max(self.last_pos, int(pos))

    def _maybe_rotate_segment(self) -> None:
        if self.segment_size <= 0:
            return
        if (
            self.inherited_segment_lines >= self.segment_size
            or self.mendelian_bad_segment_lines >= self.segment_size
        ):
            self._finish_segment()

    def _finish_segment(self) -> None:
        if self._inherited is not None:
            self._inherited.close()
            self._inherited = None
        if self._mendelian_bad is not None:
            self._mendelian_bad.close()
            self._mendelian_bad = None

        if self.segment_inherited_per_variant and self.segment_size > 0:
            write_json_atomic(
                self.output_dir / f"inherited_per_variant_seg{self.segment_index:05d}.json",
                self.segment_inherited_per_variant,
            )
            self.segment_inherited_per_variant.clear()

        self._write_cumulative_stats()
        save_checkpoint(
            self.output_dir,
            Checkpoint(
                chrom=self.last_chrom,
                last_pos=self.last_pos,
                segment_index=self.segment_index,
                cumulative=self.cumulative,
                completed=False,
            ),
        )

        self.segment_index += 1
        self._open_segment_writers()

    def close(self) -> None:
        if self._inherited is not None:
            self._inherited.close()
            self._inherited = None
        if self._mendelian_bad is not None:
            self._mendelian_bad.close()
            self._mendelian_bad = None

        if self.segment_inherited_per_variant and self.segment_size > 0:
            write_json_atomic(
                self.output_dir / f"inherited_per_variant_seg{self.segment_index:05d}.json",
                self.segment_inherited_per_variant,
            )
            self.segment_inherited_per_variant.clear()

    def finalize(self, *, completed: bool = True) -> None:
        self._write_cumulative_stats()
        if self.segment_size > 0 or completed:
            save_checkpoint(
                self.output_dir,
                Checkpoint(
                    chrom=self.last_chrom,
                    last_pos=self.last_pos,
                    segment_index=self.segment_index,
                    cumulative=self.cumulative,
                    completed=completed,
                ),
            )
        self._merge_inherited_per_variant_files()
        self.save_summary_files()

    def _write_cumulative_stats(self) -> None:
        write_json_atomic(
            self.output_dir / STATS_CUMULATIVE_FILENAME,
            {
                **self.cumulative.to_dict(),
                "last_chrom": self.last_chrom,
                "last_pos": self.last_pos,
                "segment_index": self.segment_index,
            },
        )

    def _merge_inherited_per_variant_files(self) -> None:
        merged: dict[str, int] = {}
        for path in sorted(self.output_dir.glob("inherited_per_variant_seg*.json")):
            with path.open(encoding="utf-8") as handle:
                merged.update(json.load(handle))
        if self.segment_inherited_per_variant:
            merged.update(self.segment_inherited_per_variant)
        write_json_atomic(self.output_dir / "inherited_per_variant.json", merged)

    def save_summary_files(self) -> None:
        write_json_atomic(
            self.output_dir / "inherited_per_person.json",
            self.cumulative.inherited_per_person,
        )
        write_json_atomic(
            self.output_dir / "mendelian_bad_per_gt.json",
            self.cumulative.mendelian_bad_per_gt,
        )
        write_json_atomic(
            self.output_dir / "stats.json",
            {
                "variants_seen": self.cumulative.variants_seen,
                "alleles_tested": self.cumulative.alleles_tested,
                "inherited_entries": self.cumulative.inherited_entries,
                "inherited_variants": self.cumulative.inherited_variants,
                "mendelian_bad_entries": self.cumulative.mendelian_bad_entries,
                "mendelian_bad_variants": self.cumulative.mendelian_bad_variants,
            },
        )

    @property
    def inherited_entries(self) -> int:
        return self.cumulative.inherited_entries

    @property
    def mendelian_bad_entries(self) -> int:
        return self.cumulative.mendelian_bad_entries

    @property
    def inherited_variants(self) -> int:
        return self.cumulative.inherited_variants

    @property
    def mendelian_bad_variants(self) -> int:
        return self.cumulative.mendelian_bad_variants


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


def glob_result_tsvs(output_dir: Path, prefix: str) -> list[Path]:
    segmented = sorted(output_dir.glob(f"{prefix}_*.tsv"))
    if segmented:
        return segmented
    single = output_dir / f"{prefix}.tsv"
    return [single] if single.is_file() else []
