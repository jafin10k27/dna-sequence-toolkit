"""
Nucleotide Composition and GC Analysis Module
==============================================
Calculates sequence lengths, nucleotide counts, percentages, and GC-content
with explicit reconciliation and zero-division safety.

Biological Significance of GC-Content
--------------------------------------
DNA is a double helix held together by base pairing:
  - Adenine (A) pairs with Thymine (T) via 2 hydrogen bonds.
  - Guanine (G) pairs with Cytosine (C) via 3 hydrogen bonds.

Because G-C pairs have three hydrogen bonds, DNA regions with high GC-content:
  1. Have higher melting temperatures (Tm), requiring higher denaturation temperatures
     in Polymerase Chain Reaction (PCR) assays.
  2. Exhibit greater mechanical and thermal stability.
  3. Correlate with gene density and functional genomic regions in many eukaryotes
     (e.g., CpG islands in mammalian promoter regions).

Reconciliation and Percentage Calculations
------------------------------------------
In Strict mode:
    Total Length = Count(A) + Count(T) + Count(G) + Count(C)
    GC % = ((Count(G) + Count(C)) / Total Length) * 100.0

In Ambiguous mode:
    Total Length = Count(Canonical) + Count(Degenerate Bases)
    Overall GC % = ((Count(G) + Count(C)) / Total Length) * 100.0
    Unambiguous GC % = ((Count(G) + Count(C)) / (Count(A) + Count(T) + Count(G) + Count(C))) * 100.0
    (This ensures ambiguous bases never artificially dilute GC calculations without disclosure).
"""

from collections import Counter
from dataclasses import dataclass, field
from typing import Dict, Union

from dna_toolkit.parser import DnaRecord
from dna_toolkit.validation import (
    CANONICAL_BASES,
    EmptySequenceError,
    ValidationMode,
    validate_sequence,
)


@dataclass(frozen=True)
class CompositionResult:
    """
    Structured outcome of nucleotide composition and GC-content analysis.

    Attributes:
        record_id: Identifier of the analyzed record.
        length: Total sequence length (nucleotide count).
        counts: Dictionary of individual nucleotide frequencies (e.g. {'A': 3, 'T': 2, ...}).
        percentages: Percentage of each base relative to total length (0.0 to 100.0).
        gc_count: Total number of G and C nucleotides.
        gc_percent: Percentage of G and C bases relative to total sequence length.
        ambiguous_count: Number of non-canonical (IUPAC degenerate) bases.
        unambiguous_gc_percent: GC percentage computed only over canonical (A, T, G, C) bases.
    """
    record_id: str
    length: int
    counts: Dict[str, int]
    percentages: Dict[str, float]
    gc_count: int
    gc_percent: float
    ambiguous_count: int = 0
    unambiguous_gc_percent: float = 0.0


def calculate_gc_content(sequence: str) -> float:
    """
    Calculate the GC-content percentage for a validated DNA sequence.

    Args:
        sequence: DNA sequence string (case-insensitive).

    Returns:
        GC percentage as a float between 0.0 and 100.0.

    Raises:
        EmptySequenceError: If the sequence is empty.
    """
    cleaned = "".join(sequence.split()).upper()
    if not cleaned:
        raise EmptySequenceError("Cannot calculate GC-content of an empty sequence.")

    gc_count = cleaned.count("G") + cleaned.count("C")
    return (gc_count / len(cleaned)) * 100.0


def analyze_composition(
    record: Union[DnaRecord, str],
    mode: Union[ValidationMode, str] = ValidationMode.STRICT,
    record_id: str = "unnamed_record",
) -> CompositionResult:
    """
    Perform comprehensive nucleotide frequency and GC-content analysis.

    Args:
        record: A DnaRecord instance or a raw DNA sequence string.
        mode: Validation policy (STRICT or AMBIGUOUS).
        record_id: Identifier used if record is passed as a string.

    Returns:
        CompositionResult containing lengths, counts, percentages, and GC metrics.

    Raises:
        EmptySequenceError: If the sequence has length 0.
        InvalidBaseError: If the sequence contains characters disallowed under mode.
        ValueError: If composition totals fail reconciliation with sequence length.
    """
    if isinstance(record, DnaRecord):
        rec_id = record.id
        raw_seq = record.sequence
    else:
        rec_id = record_id
        raw_seq = record

    # Validate and normalize sequence first
    valid_seq = validate_sequence(raw_seq, mode=mode, record_id=rec_id)
    seq_len = len(valid_seq)

    if seq_len == 0:
        raise EmptySequenceError(f"Sequence for record '{rec_id}' is empty.")

    counter = Counter(valid_seq)

    # Canonical base counts
    count_a = counter.get("A", 0)
    count_t = counter.get("T", 0)
    count_g = counter.get("G", 0)
    count_c = counter.get("C", 0)

    gc_count = count_g + count_c
    gc_percent = (gc_count / seq_len) * 100.0

    # Collect all counts sorted canonically then alphabetically
    ordered_counts: Dict[str, int] = {
        "A": count_a,
        "T": count_t,
        "G": count_g,
        "C": count_c,
    }
    
    ambiguous_count = 0
    for base, count in counter.items():
        if base not in CANONICAL_BASES:
            ordered_counts[base] = count
            ambiguous_count += count

    # Reconciliation integrity check
    total_reconciled = sum(ordered_counts.values())
    if total_reconciled != seq_len:
        raise ValueError(
            f"Composition reconciliation error for '{rec_id}': "
            f"sum of nucleotide counts ({total_reconciled}) does not match sequence length ({seq_len})."
        )

    # Compute percentages
    percentages: Dict[str, float] = {
        base: round((cnt / seq_len) * 100.0, 4)
        for base, cnt in ordered_counts.items()
    }

    # Unambiguous GC percentage
    canonical_total = count_a + count_t + count_g + count_c
    unambiguous_gc_percent = (
        (gc_count / canonical_total * 100.0) if canonical_total > 0 else 0.0
    )

    return CompositionResult(
        record_id=rec_id,
        length=seq_len,
        counts=ordered_counts,
        percentages=percentages,
        gc_count=gc_count,
        gc_percent=round(gc_percent, 4),
        ambiguous_count=ambiguous_count,
        unambiguous_gc_percent=round(unambiguous_gc_percent, 4),
    )
