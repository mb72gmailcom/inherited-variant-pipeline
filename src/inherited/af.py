from __future__ import annotations

import json
from pathlib import Path

from inherited.constants import DEFAULT_AF_THRESHOLD


def load_af_json(path: Path) -> dict[str, float]:
    """Load gnomAD allele frequencies from JSON.

    Expected format::

        {"22:12345:A:G": 0.0001, "var_key": {"AF": 0.001, "AF_EUR": 0.0005}}

    For object values, AF_EUR is preferred when present, otherwise AF.
  Missing keys default to 0 at lookup time (variant kept).
    """
    with path.open(encoding="utf-8") as handle:
        raw = json.load(handle)

    if not isinstance(raw, dict):
        raise ValueError(f"AF JSON must be an object at top level: {path}")

    table: dict[str, float] = {}
    for key, value in raw.items():
        table[str(key)] = _af_value(value)
    return table


def _af_value(value: object) -> float:
    if value is None:
        return 0.0
    if isinstance(value, (int, float)):
        return float(value)
    if isinstance(value, dict):
        if value.get("AF_EUR") not in (None, "", "."):
            return float(value["AF_EUR"])
        if value.get("AF") not in (None, "", "."):
            return float(value["AF"])
        return 0.0
    raise ValueError(f"Invalid AF JSON value: {value!r}")


def is_rare(af_table: dict[str, float], key: str, threshold: float = DEFAULT_AF_THRESHOLD) -> bool:
    return af_table.get(key, 0.0) <= threshold
