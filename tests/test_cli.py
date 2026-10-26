"""
Tests for Command-Line Interface (CLI)
======================================
Validates arguments, operation flags, output formats, and process exit codes:
  0: Success
  1: Data validation or parsing error
  2: Command line argument error
  3: I/O or filesystem error
"""

from pathlib import Path
import pytest

from dna_toolkit.cli import (
    EXIT_ARG_ERROR,
    EXIT_DATA_ERROR,
    EXIT_IO_ERROR,
    EXIT_SUCCESS,
    main,
)


def test_cli_help(capsys):
    """Test that --help prints usage and exits with code 0."""
    with pytest.raises(SystemExit) as exc_info:
        main(["--help"])
    assert exc_info.value.code == 0
    captured = capsys.readouterr()
    assert "DNA Sequence Analysis Toolkit" in captured.out


def test_cli_version(capsys):
    """Test that --version prints version and exits with code 0."""
    with pytest.raises(SystemExit) as exc_info:
        main(["--version"])
    assert exc_info.value.code == 0
    captured = capsys.readouterr()
    assert "dna-toolkit" in captured.out


def test_cli_missing_input_argument(capsys):
    """Test running CLI without specifying input FASTA file returns EXIT_ARG_ERROR."""
    exit_code = main([])
    assert exit_code == EXIT_ARG_ERROR
    captured = capsys.readouterr()
    assert "Missing input FASTA file" in captured.err


def test_cli_nonexistent_file(capsys, tmp_path: Path):
    """Test running CLI on non-existent file returns EXIT_IO_ERROR (code 3)."""
    fake_path = str(tmp_path / "not_there.fasta")
    exit_code = main([fake_path])
    assert exit_code == EXIT_IO_ERROR
    captured = capsys.readouterr()
    assert "File Error" in captured.err


def test_cli_empty_fasta_file(empty_fasta_file: Path, capsys):
    """Test running CLI on empty FASTA file returns EXIT_DATA_ERROR (code 1)."""
    exit_code = main([str(empty_fasta_file)])
    assert exit_code == EXIT_DATA_ERROR
    captured = capsys.readouterr()
    assert "FASTA Format Error" in captured.err


def test_cli_invalid_dna_strict_mode(capsys, tmp_path: Path):
    """Test running CLI with invalid sequence in strict mode returns EXIT_DATA_ERROR (code 1)."""
    corrupt_file = tmp_path / "corrupt.fasta"
    corrupt_file.write_text(">bad_seq\nATGCZATGC\n", encoding="utf-8")

    exit_code = main([str(corrupt_file)])
    assert exit_code == EXIT_DATA_ERROR
    captured = capsys.readouterr()
    assert "Validation Error" in captured.err
    assert "'Z'" in captured.err


def test_cli_duplicate_ids_without_and_with_flag(duplicate_ids_fasta_file: Path, capsys):
    """Test duplicate ID handling without and with --allow-duplicates flag."""
    # Without flag -> error
    exit_code = main([str(duplicate_ids_fasta_file)])
    assert exit_code == EXIT_DATA_ERROR
    captured = capsys.readouterr()
    assert "Duplicate sequence identifier" in captured.err

    # With flag -> success
    exit_code_ok = main([str(duplicate_ids_fasta_file), "--allow-duplicates"])
    assert exit_code_ok == EXIT_SUCCESS


def test_cli_ambiguous_mode_flag(ambiguous_fasta_file: Path, capsys):
    """Test that ambiguous FASTA succeeds with --mode ambiguous."""
    # Strict fails
    exit_strict = main([str(ambiguous_fasta_file), "--mode", "strict"])
    assert exit_strict == EXIT_DATA_ERROR

    # Ambiguous succeeds
    exit_ambig = main([str(ambiguous_fasta_file), "--mode", "ambiguous"])
    assert exit_ambig == EXIT_SUCCESS
    captured = capsys.readouterr()
    assert "Ambiguous Bases" in captured.out


def test_cli_selective_operations(sample_fasta_file: Path, capsys):
    """Test running only composition analysis."""
    exit_code = main([str(sample_fasta_file), "--composition"])
    assert exit_code == EXIT_SUCCESS
    captured = capsys.readouterr()
    assert "GC-Content" in captured.out
    # Translation should not have run because only --composition was selected
    assert "Translation:" not in captured.out


def test_cli_json_export_to_file(sample_fasta_file: Path, tmp_path: Path, capsys):
    """Test exporting JSON to file."""
    out_json = tmp_path / "out.json"
    exit_code = main([str(sample_fasta_file), "--format", "json", "-o", str(out_json)])
    assert exit_code == EXIT_SUCCESS
    assert out_json.is_file()
    captured = capsys.readouterr()
    assert "Successfully exported JSON" in captured.out


def test_cli_csv_export_to_file(sample_fasta_file: Path, tmp_path: Path, capsys):
    """Test exporting CSV to file."""
    out_csv = tmp_path / "out.csv"
    exit_code = main([str(sample_fasta_file), "--format", "csv", "-o", str(out_csv)])
    assert exit_code == EXIT_SUCCESS
    assert out_csv.is_file()
    captured = capsys.readouterr()
    assert "Successfully exported CSV" in captured.out


def test_cli_translation_options(sample_fasta_file: Path, capsys):
    """Test passing custom reading frame and genetic code table."""
    exit_code = main([
        str(sample_fasta_file),
        "--translate",
        "--frame", "2",
        "--table", "2",
    ])
    assert exit_code == EXIT_SUCCESS
    captured = capsys.readouterr()
    assert "Frame: +2" in captured.out
    assert "Code Table: 2" in captured.out
