"""
DNA Sequence Analysis Toolkit
==============================
A modular, educational bioinformatics toolkit built for parsing, validating,
analyzing, and translating DNA sequences from FASTA files.

Designed for students and engineers learning computational biology and Python.
"""

from dna_toolkit.parser import DnaRecord, parse_fasta
from dna_toolkit.validation import (
    ValidationMode,
    validate_sequence,
    validate_records,
    ValidationError,
    EmptySequenceError,
    InvalidBaseError,
    DuplicateRecordError,
)
from dna_toolkit.analysis import (
    CompositionResult,
    analyze_composition,
    calculate_gc_content,
)
from dna_toolkit.operations import (
    complement,
    reverse_complement,
    translate_dna,
    TranslationResult,
)

__version__ = "0.1.0"
__all__ = [
    "DnaRecord",
    "parse_fasta",
    "ValidationMode",
    "validate_sequence",
    "validate_records",
    "ValidationError",
    "EmptySequenceError",
    "InvalidBaseError",
    "DuplicateRecordError",
    "CompositionResult",
    "analyze_composition",
    "calculate_gc_content",
    "complement",
    "reverse_complement",
    "translate_dna",
    "TranslationResult",
    "__version__",
]
