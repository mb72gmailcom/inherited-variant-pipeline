from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class FamilyRelations:
    """Family relation tables built from the family TSV file."""

    trio: dict[str, tuple[str, str]] = field(default_factory=dict)
    trio_cl: dict[str, tuple[str, str]] = field(default_factory=dict)
    trio_all: dict[str, tuple[str, str]] = field(default_factory=dict)
    counts: dict[str, int] = field(default_factory=dict)
    family_size: dict[str, int] = field(default_factory=dict)
    trios_ids: list[list[str]] = field(default_factory=list)


def load_family_relations(path: Path) -> FamilyRelations:
    """Load family relations from a tab-separated file.

    Expected columns (no header name requirement beyond first cell ``spid``)::

        spid    family_id    father_id    mother_id

    ``dTrioCl`` maps child sample id -> (mother_id, father_id).
    """
    relations = FamilyRelations()

    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row or row[0] == "spid":
                continue

            spid, family_id, father_id, mother_id = row[0], row[1], row[2], row[3]
            relations.family_size[family_id] = relations.family_size.get(family_id, 0) + 1
            relations.counts[spid] = relations.counts.get(spid, 0) + 1
            relations.trio_all[spid] = (father_id, mother_id)

            if father_id != "0" or mother_id != "0":
                relations.trio[spid] = (mother_id, father_id)

            if father_id != "0" and mother_id != "0":
                relations.trio_cl[spid] = (mother_id, father_id)
                relations.trios_ids.append([spid, father_id, mother_id])

    return relations


def build_trio_indices(
    sample_header: list[str],
    trio_cl: dict[str, tuple[str, str]],
) -> tuple[dict[int, tuple[int, int]], list[tuple[int, int, int]]]:
    """Map VCF column indices for complete child-mother-father trios."""
    pid_to_idx = {pid: i for i, pid in enumerate(sample_header)}
    trio_ind: dict[int, tuple[int, int]] = {}
    trios_ind: list[tuple[int, int, int]] = []

    for child_idx, pid in enumerate(sample_header):
        if pid not in trio_cl:
            continue
        mother_id, father_id = trio_cl[pid]
        mother_idx = pid_to_idx.get(mother_id)
        father_idx = pid_to_idx.get(father_id)
        if mother_idx is None or father_idx is None:
            continue
        trio_ind[child_idx] = (mother_idx, father_idx)
        trios_ind.append((child_idx, mother_idx, father_idx))

    return trio_ind, trios_ind
