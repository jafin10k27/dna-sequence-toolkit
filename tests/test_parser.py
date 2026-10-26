"""
Tests for FASTA Parser
======================
Validates FASTA reading, multiline wrapping, empty files, missing files,
malformed formats, and character preservation.
"""

from pathlib import Path
import pytest

from dna_toolkit.parser import (
    DnaRecord,
    EmptyFastaError,
    FastaFileNotFoundError,
    MalformedFastaError,
    parse_fasta,
)


def test_parse_single_record(tmp_path: Path):
    """Test parsing a single FASTA record."""
    file = tmp_path / "single.fasta"
    file.write_text(">seq_01 Homo sapiens test\nACGTACGT\n", encoding="utf-8")

    records = parse_fasta(file)
    assert len(records) == 1
    assert records[0].id == "seq_01"
    assert records[0].description == "seq_01 Homo sapiens test"
    assert records[0].sequence == "ACGTACGT"


def test_parse_multi_record(sample_fasta_file: Path):
    """Test parsing a file containing multiple records."""
    records = parse_fasta(sample_fasta_file)
    assert len(records) == 3
    assert records[0].id == "seq1"
    assert records[0].sequence == "ATGC"
    assert records[1].id == "seq2"
    assert records[1].sequence == "GATTACAGATTACA"
    assert records[2].id == "seq3"
    assert records[2].sequence == "ATGGAA"


def test_parse_multiline_sequence(multiline_fasta_file: Path):
    """Test that multiline sequences are concatenated without loss or spurious characters."""
    records = parse_fasta(multiline_fasta_file)
    assert len(records) == 1
    expected_sequence = "ATGCATGCGGCCGGCCTTTTAAAA"
    assert records[0].sequence == expected_sequence
    assert len(records[0].sequence) == 24


def test_parse_missing_file(tmp_path: Path):
    """Test that attempting to parse a non-existent file raises FastaFileNotFoundError."""
    missing_path = tmp_path / "does_not_exist.fasta"
    with pytest.raises(FastaFileNotFoundError) as exc_info:
        parse_fasta(missing_path)
    assert "not found" in str(exc_info.value).lower()


def test_parse_empty_file(empty_fasta_file: Path):
    """Test that an empty 0-byte file raises EmptyFastaError."""
    with pytest.raises(EmptyFastaError) as exc_info:
        parse_fasta(empty_fasta_file)
    assert "empty" in str(exc_info.value).lower()


def test_parse_whitespace_only_file(whitespace_fasta_file: Path):
    """Test that a file with only spaces and newlines raises EmptyFastaError."""
    with pytest.raises(EmptyFastaError) as exc_info:
        parse_fasta(whitespace_fasta_file)
    assert "whitespace" in str(exc_info.value).lower()


def test_parse_non_fasta_file(non_fasta_file: Path):
    """Test that a text file lacking FASTA '>' headers raises MalformedFastaError."""
    with pytest.raises(MalformedFastaError) as exc_info:
        parse_fasta(non_fasta_file)
    assert "header starting with '>'" in str(exc_info.value).lower()


def test_parser_does_not_silently_discard_invalid_characters(tmp_path: Path):
    """Test that parser does not silently strip or alter unexpected sequence characters."""
    file = tmp_path / "corrupt.fasta"
    file.write_text(">bad_seq Corrupt entry\nATG123*XYZ\n", encoding="utf-8")

    records = parse_fasta(file)
    assert len(records) == 1
    # Parser must preserve raw characters for validation to inspect
    assert records[0].sequence == "ATG123*XYZ"
