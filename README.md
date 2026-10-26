# DNA Sequence Toolkit

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![Tests](https://img.shields.io/badge/tests-50%20passed-brightgreen.svg)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Code Style: Black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

A modular command-line tool and Python library for parsing, validating, analyzing, and translating DNA sequences from FASTA files.

Built on Biopython, `dna-toolkit` provides strict and degenerate (IUPAC) sequence validation, GC-content and nucleotide frequency metrics, strand complementation, 6-frame translation with partial codon handling, and multi-format exports (Terminal, JSON, and CSV).

---

## Features

- **FASTA Parsing**: Handles single and multi-record files with multiline sequence wrapping via Biopython `SeqIO`.
- **Validation Policies**:
  - `strict`: Enforces canonical DNA nucleotides (`A`, `T`, `G`, `C`) with case normalization.
  - `ambiguous`: Permits IUPAC degenerate codes (`R`, `Y`, `S`, `W`, `K`, `M`, `B`, `D`, `H`, `V`, `N`).
  - Strict duplicate header rejection with optional override (`--allow-duplicates`).
- **Nucleotide & GC Analysis**:
  - Calculates length, nucleotide counts, percentages, and total GC%.
  - Reconciles composition totals against total length.
  - Computes unambiguous GC-content when ambiguous bases are present.
- **Sequence Operations**:
  - Direct strand complement ($3' \to 5'$) and reverse complement ($5' \to 3'$).
  - 6-frame translation (+1, +2, +3 forward; -1, -2, -3 reverse strand).
  - NCBI genetic code table support (Standard, Vertebrate Mitochondrial, Bacterial, etc.).
  - Automatic trimming and accounting for partial trailing codons.
- **Reporting & Exports**:
  - Human-readable terminal summaries.
  - Structured JSON export with execution metadata.
  - RFC-4180-compliant CSV export with normalized sequence fields.

---

## Installation

### From Source

Clone the repository and install in a virtual environment:

```bash
git clone https://github.com/jafin10k27/dna-sequence-toolkit.git
cd dna-sequence-toolkit

# Create and activate virtual environment
python3 -m venv .venv
source .venv/bin/activate       # Bash / Zsh
# source .venv/bin/activate.fish # Fish shell

# Install package and test dependencies
pip install -e ".[dev]"
```

---

## Quickstart

### Command-Line Interface

Run a full analysis on a FASTA file:
```bash
dna-toolkit data/sample.fasta
```

#### Selective Operations
```bash
# Nucleotide frequency and GC-content only
dna-toolkit data/sample.fasta --composition

# Generate reverse complement only
dna-toolkit data/sample.fasta --reverse-complement

# Translate in reading frame 2 using Vertebrate Mitochondrial code (NCBI table 2)
dna-toolkit data/sample.fasta --translate --frame 2 --table 2

# Halt translation at the first stop codon
dna-toolkit data/sample.fasta --translate --to-stop

# Validate files containing IUPAC ambiguous bases (e.g. N, R, Y)
dna-toolkit data/ambiguous.fasta --mode ambiguous
```

#### Exporting Results
```bash
# Export full analysis to JSON
dna-toolkit data/sample.fasta --format json -o output/report.json

# Export to CSV
dna-toolkit data/sample.fasta --format csv -o output/report.csv
```

---

### Python API

`dna-toolkit` can also be imported directly into Python scripts and pipelines:

```python
from dna_toolkit import parse_fasta, analyze_composition, reverse_complement, translate_dna

# 1. Parse records from FASTA
records = parse_fasta("data/sample.fasta")

for record in records:
    print(f"Record: {record.id}")

    # 2. Composition and GC-content
    comp = analyze_composition(record.sequence)
    print(f"  Length: {comp.length} bp")
    print(f"  GC: {comp.gc_percent:.2f}% (G+C: {comp.gc_count})")
    print(f"  Counts: {comp.counts}")

    # 3. Reverse complement (5' to 3')
    rev = reverse_complement(record.sequence)
    print(f"  RevComp: {rev[:20]}...")

    # 4. Translation
    trans = translate_dna(record.sequence, frame=1, table=1)
    print(f"  Protein: {trans.protein_sequence[:20]}... ({trans.codon_count} aa)")
```

---

## CLI Reference

```
usage: dna-toolkit [-h] [-i INPUT_FILE] [-v] [--all] [--composition]
                   [--reverse-complement] [--complement] [--translate]
                   [-m {strict,ambiguous}] [--allow-duplicates]
                   [-f {1,2,3,-1,-2,-3}] [-t TABLE] [--to-stop]
                   [--format {text,json,csv}] [-o OUTPUT]
                   [fasta_file]
```

### Options

| Flag | Description | Default |
| :--- | :--- | :--- |
| `fasta_file` / `-i` | Path to input FASTA file | Positional / required |
| `--all` | Run all analysis and translation operations | Active if no specific op given |
| `--composition` | Calculate sequence length, frequencies, and GC% | Off |
| `--reverse-complement`| Generate $5' \to 3'$ reverse complement | Off |
| `--complement` | Generate direct $3' \to 5'$ strand complement | Off |
| `--translate` | Translate DNA sequence to protein | Off |
| `-m`, `--mode` | Validation mode: `strict` or `ambiguous` | `strict` |
| `--allow-duplicates` | Permit duplicate record IDs in input | `False` |
| `-f`, `--frame` | Translation reading frame (`1, 2, 3, -1, -2, -3`) | `1` |
| `-t`, `--table` | NCBI genetic code table ID or name | `1` (Standard) |
| `--to-stop` | Halt translation at first stop codon | `False` |
| `--format` | Output format: `text`, `json`, or `csv` | `text` |
| `-o`, `--output` | Destination file path (prints to stdout if omitted) | `None` |

### Exit Codes

| Code | Meaning |
| :--- | :--- |
| `0` | Success (`EXIT_SUCCESS`) |
| `1` | Data validation or parsing error (`EXIT_DATA_ERROR`) |
| `2` | Invalid CLI argument or option syntax (`EXIT_ARG_ERROR`) |
| `3` | File not found or I/O permission error (`EXIT_IO_ERROR`) |

---

## Project Structure

```
dna-sequence-toolkit/
├── src/
│   └── dna_toolkit/
│       ├── __init__.py       # Public package API
│       ├── parser.py         # FASTA reader and DnaRecord model
│       ├── validation.py     # Strict and ambiguous validation
│       ├── analysis.py       # Composition and GC calculation
│       ├── operations.py     # Complementation and 6-frame translation
│       ├── reporting.py      # Terminal, JSON, and CSV exporters
│       └── cli.py            # CLI entry point and argument parsing
├── tests/                    # 50 unit and integration tests
├── data/                     # Sample datasets (valid, multiline, ambiguous, edge cases)
├── docs/
│   ├── architecture.md       # Pipeline data flow and design choices
│   └── interview_guide.md    # Technical interview questions and domain concepts
├── pyproject.toml            # Build configuration and dependency specification
└── main.py                   # Direct script entry point
```

---

## Testing

Run tests with `pytest`:

```bash
pytest -v
```

The test suite covers:
- Single-record, multi-record, and multiline FASTA parsing.
- Empty files, whitespace-only files, and missing files.
- Strict and ambiguous IUPAC nucleotide validation.
- Lowercase sequence normalization.
- Incomplete codon trimming and stop codon behaviors.
- Reading frame selection across forward (+1, +2, +3) and reverse (-1, -2, -3) strands.
- JSON and CSV output integrity.
- CLI argument parsing, flags, and exit codes.

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
