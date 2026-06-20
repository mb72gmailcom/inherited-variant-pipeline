from __future__ import annotations

import resource
import sys
from pathlib import Path

from inherited.constants import DEFAULT_MEMORY_BLOCK


def get_memory_usage_mb() -> float:
    """Return resident memory usage in MiB."""
    proc_status = Path("/proc/self/status")
    if proc_status.is_file():
        for line in proc_status.read_text(encoding="utf-8").splitlines():
            if line.startswith("VmRSS:"):
                return int(line.split()[1]) / 1024.0

    usage = resource.getrusage(resource.RUSAGE_SELF)
    if sys.platform == "darwin":
        return usage.ru_maxrss / (1024 * 1024)
    return usage.ru_maxrss / 1024.0


def log_memory_if_due(
    variants_seen: int,
    *,
    debug: bool,
    memory_block: int = DEFAULT_MEMORY_BLOCK,
) -> None:
    """Print memory usage every ``memory_block`` processed variants when debug is on."""
    if not debug or memory_block <= 0:
        return
    if variants_seen > 0 and variants_seen % memory_block == 0:
        memory_mb = get_memory_usage_mb()
        print(
            f"[debug] variants={variants_seen:,} memory={memory_mb:.1f} MiB",
            flush=True,
        )
