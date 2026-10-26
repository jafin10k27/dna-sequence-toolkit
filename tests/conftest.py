"""
Pytest Fixtures and Test Helpers
================================
Provides reusable fixtures for temporary FASTA files, mock inputs, and shared test data.
"""

from pathlib import Path
import pytest


@pytest.fixture
def sample_fasta_file(tmp_path: Path) -> Path:
    """Fixture creating a multi-record FASTA file with canonical sequences."""
    content = (
        ">seq1 First record description\n"
        "ATGC\n"
        ">seq2 Second multiline record\n"
        "GATTACA\n"
        "GATTACA\n"
        ">seq3 Coding test\n"
        "ATGGAA\n"
    )
    fasta_path = tmp_path / "test_sample.fasta"
    fasta_path.write_text(content, encoding="utf-8")
    return fasta_path


@pytest.fixture
def multiline_fasta_file(tmp_path: Path) -> Path:
    """Fixture creating a single record wrapped across multiple lines."""
    content = (
        ">multiline_record A long multiline sequence\n"
        "ATGCATGC\n"
        "GGCCGGCC\n"
        "TTTTAAAA\n"
    )
    fasta_path = tmp_path / "test_multiline.fasta"
    fasta_path.write_text(content, encoding="utf-8")
    return fasta_path


@pytest.fixture
def empty_fasta_file(tmp_path: Path) -> Path:
    """Fixture creating an empty 0-byte file."""
    fasta_path = tmp_path / "test_empty.fasta"
    fasta_path.write_text("", encoding="utf-8")
    return fasta_path


@pytest.fixture
def whitespace_fasta_file(tmp_path: Path) -> Path:
    """Fixture creating a file containing only whitespace."""
    fasta_path = tmp_path / "test_whitespace.fasta"
    fasta_path.write_text("   \n\t\n  \n", encoding="utf-8")
    return fasta_path


@pytest.fixture
def non_fasta_file(tmp_path: Path) -> Path:
    """Fixture creating a plain text file without FASTA headers."""
    txt_path = tmp_path / "test_plain.txt"
    txt_path.write_text("This is a plain text file without FASTA formatting.\n", encoding="utf-8")
    return txt_path


@pytest.fixture
def ambiguous_fasta_file(tmp_path: Path) -> Path:
    """Fixture creating a FASTA file with IUPAC degenerate codes."""
    content = (
        ">ambig1 Ambiguous sequence\n"
        "ATGCRAMYKSWN\n"
    )
    fasta_path = tmp_path / "test_ambig.fasta"
    fasta_path.write_text(content, encoding="utf-8")
    return fasta_path


@pytest.fixture
def duplicate_ids_fasta_file(tmp_path: Path) -> Path:
    """Fixture creating a FASTA file with duplicate record identifiers."""
    content = (
        ">dup_id First\n"
        "ATGC\n"
        ">dup_id Second\n"
        "CGTA\n"
    )
    fasta_path = tmp_path / "test_dups.fasta"
    fasta_path.write_text(content, encoding="utf-8")
    return fasta_path
