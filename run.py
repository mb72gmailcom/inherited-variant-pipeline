#!/usr/bin/env python3
"""Run the inherited CLI without installing the package.

Usage:
    python run.py analyze --vcf file.vcf.gz --af-json af.json -o results/
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from inherited.cli import main

if __name__ == "__main__":
    main()
