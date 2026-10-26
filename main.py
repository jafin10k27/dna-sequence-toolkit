#!/usr/bin/env python3
"""
Top-level entry script for running the DNA Sequence Analysis Toolkit.

Usage:
    python main.py --help
    python main.py data/sample.fasta
"""

import sys
from pathlib import Path

# Ensure src layout is available in Python path even without pip install -e .
PROJECT_ROOT = Path(__file__).resolve().parent
SRC_PATH = PROJECT_ROOT / "src"
if str(SRC_PATH) not in sys.path:
    sys.path.insert(0, str(SRC_PATH))

from dna_toolkit.cli import main

if __name__ == "__main__":
    sys.exit(main())
