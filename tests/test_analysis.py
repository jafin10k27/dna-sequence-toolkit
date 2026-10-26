"""
Tests for Nucleotide Composition and GC-Content Analysis
========================================================
Validates known biological benchmarks, composition reconciliation, nucleotide
frequency counts, percentages, and ambiguous mode metrics.
"""

import pytest

from dna_toolkit.analysis import analyze_composition, calculate_gc_content
from dna_toolkit.parser import DnaRecord
from dna_toolkit.validation import EmptySequenceError, ValidationMode


def test_known_case_atgc():
    """Verify known case: ATGC has length 4 and GC content 50%."""
    result = analyze_composition("ATGC")
    assert result.length == 4
    assert result.gc_count == 2
    assert result.gc_percent == 50.0
    assert result.counts == {"A": 1, "T": 1, "G": 1, "C": 1}
    assert result.percentages["A"] == 25.0
    assert result.percentages["T"] == 25.0
    assert result.percentages["G"] == 25.0
    assert result.percentages["C"] == 25.0


def test_known_case_ggcc():
    """Verify known case: GGCC has GC content 100%."""
    result = analyze_composition("GGCC")
    assert result.length == 4
    assert result.gc_count == 4
    assert result.gc_percent == 100.0
    assert result.counts["G"] == 2
    assert result.counts["C"] == 2
    assert result.counts["A"] == 0
    assert result.counts["T"] == 0


def test_known_case_atat():
    """Verify known case: ATAT has GC content 0%."""
    result = analyze_composition("ATAT")
    assert result.length == 4
    assert result.gc_count == 0
    assert result.gc_percent == 0.0
    assert result.counts["A"] == 2
    assert result.counts["T"] == 2
    assert result.counts["G"] == 0
    assert result.counts["C"] == 0


def test_known_case_gattaca():
    """Verify known case: GATTACA has counts A=3, T=2, G=1, C=1."""
    result = analyze_composition("GATTACA")
    assert result.length == 7
    assert result.counts["A"] == 3
    assert result.counts["T"] == 2
    assert result.counts["G"] == 1
    assert result.counts["C"] == 1
    assert result.gc_count == 2
    # GC % = 2 / 7 * 100 = 28.5714...
    assert pytest.approx(result.gc_percent, 0.01) == 28.57


def test_composition_totals_reconciliation():
    """Verify that composition counts sum precisely to sequence length."""
    seq = "ATGCGTACGTTAGCTAGCTAGCTAAGCTAGCTA"
    result = analyze_composition(seq)
    total_counts = sum(result.counts.values())
    assert total_counts == result.length
    assert total_counts == len(seq)


def test_analyze_composition_with_record_object():
    """Verify analyze_composition works seamlessly when given a DnaRecord."""
    record = DnaRecord(id="rec01", description="test description", sequence="ATGC")
    result = analyze_composition(record)
    assert result.record_id == "rec01"
    assert result.length == 4
    assert result.gc_percent == 50.0


def test_empty_sequence_raises_error():
    """Verify that calculating GC or composition on an empty sequence raises EmptySequenceError."""
    with pytest.raises(EmptySequenceError):
        calculate_gc_content("")

    with pytest.raises(EmptySequenceError):
        analyze_composition("")


def test_ambiguous_mode_metrics():
    """Verify nucleotide analysis in ambiguous mode calculates degenerate counts and unambiguous GC%."""
    # Sequence: A, T, G, C, N, N (length 6: canonical G+C=2, canonical total=4)
    result = analyze_composition("ATGCNN", mode=ValidationMode.AMBIGUOUS)
    assert result.length == 6
    assert result.counts["N"] == 2
    assert result.ambiguous_count == 2
    # Overall GC% = 2 / 6 * 100 = 33.3333%
    assert pytest.approx(result.gc_percent, 0.01) == 33.33
    # Unambiguous GC% = 2 / 4 * 100 = 50.0%
    assert pytest.approx(result.unambiguous_gc_percent, 0.01) == 50.0
