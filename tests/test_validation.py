"""
Tests for DNA Sequence Validation
=================================
Validates strict canonical mode, ambiguous IUPAC mode, lowercase normalization,
empty sequences, invalid symbol reporting, and duplicate record handling.
"""

import pytest

from dna_toolkit.parser import DnaRecord
from dna_toolkit.validation import (
    DuplicateRecordError,
    EmptySequenceError,
    InvalidBaseError,
    ValidationMode,
    validate_records,
    validate_sequence,
)


def test_validate_canonical_strict_valid():
    """Test standard valid canonical sequences in strict mode."""
    assert validate_sequence("ATGC") == "ATGC"
    assert validate_sequence("GATTACA") == "GATTACA"
    assert validate_sequence("GGCC") == "GGCC"


def test_validate_lowercase_normalization():
    """Test that lower-case nucleotides are normalized to upper-case."""
    assert validate_sequence("atgc") == "ATGC"
    assert validate_sequence("GattAca") == "GATTACA"
    assert validate_sequence("aTgCgGcC") == "ATGCGGCC"


def test_validate_whitespace_stripping():
    """Test that internal or surrounding whitespace is cleaned."""
    assert validate_sequence("  ATG C  \n\t") == "ATGC"


def test_validate_empty_sequence_raises_error():
    """Test that an empty sequence raises EmptySequenceError."""
    with pytest.raises(EmptySequenceError) as exc_info:
        validate_sequence("", record_id="test_empty")
    assert "empty sequence" in str(exc_info.value).lower()
    assert "test_empty" in str(exc_info.value)


def test_validate_whitespace_only_sequence_raises_error():
    """Test that a sequence of only spaces/newlines raises EmptySequenceError."""
    with pytest.raises(EmptySequenceError):
        validate_sequence("   \n\t  ")


def test_validate_strict_mode_rejects_ambiguous_bases():
    """Test that strict mode rejects IUPAC ambiguous bases like N, R, Y."""
    with pytest.raises(InvalidBaseError) as exc_info:
        validate_sequence("ATGCN", mode=ValidationMode.STRICT, record_id="rec_n")
    assert "N" in exc_info.value.invalid_chars
    assert 5 in exc_info.value.char_positions["N"]


def test_validate_ambiguous_mode_accepts_iupac():
    """Test that ambiguous mode accepts valid IUPAC degenerate codes."""
    seq = "ATGCRAMYKSWNBVDH"
    validated = validate_sequence(seq, mode=ValidationMode.AMBIGUOUS)
    assert validated == seq


def test_validate_invalid_characters_detail():
    """Test that invalid symbols (e.g. '1', 'Z', '*') are identified with positions."""
    with pytest.raises(InvalidBaseError) as exc_info:
        validate_sequence("AT1GC*Z", mode=ValidationMode.STRICT)
    assert exc_info.value.invalid_chars == {"1", "*", "Z"}
    assert exc_info.value.char_positions["1"] == [3]
    assert exc_info.value.char_positions["*"] == [6]
    assert exc_info.value.char_positions["Z"] == [7]


def test_validate_records_duplicate_identifiers_disallowed():
    """Test that duplicate record IDs raise DuplicateRecordError by default."""
    records = [
        DnaRecord(id="seq1", description="desc1", sequence="ATGC"),
        DnaRecord(id="seq1", description="desc2", sequence="CGTA"),
    ]
    with pytest.raises(DuplicateRecordError) as exc_info:
        validate_records(records, allow_duplicates=False)
    assert "Duplicate sequence identifier 'seq1'" in str(exc_info.value)


def test_validate_records_duplicate_identifiers_allowed():
    """Test that duplicate IDs succeed when allow_duplicates=True."""
    records = [
        DnaRecord(id="seq1", description="desc1", sequence="ATGC"),
        DnaRecord(id="seq1", description="desc2", sequence="CGTA"),
    ]
    validated = validate_records(records, allow_duplicates=True)
    assert len(validated) == 2
    assert validated[0].sequence == "ATGC"
    assert validated[1].sequence == "CGTA"
