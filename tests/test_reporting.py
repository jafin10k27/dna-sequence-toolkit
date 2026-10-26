"""
Tests for Reporting and Export Module
=====================================
Validates human-readable terminal reports, JSON structure, and RFC-4180 CSV compliance.
"""

import csv
import json
from pathlib import Path
import pytest

from dna_toolkit.analysis import analyze_composition
from dna_toolkit.operations import reverse_complement, translate_dna
from dna_toolkit.reporting import (
    RecordReport,
    export_csv,
    export_json,
    format_terminal_report,
)


@pytest.fixture
def sample_reports():
    """Build a pair of RecordReport fixtures for testing exporters."""
    seq1 = "ATGC"
    r1 = RecordReport(
        record_id="seq1",
        description="seq1 First test",
        sequence=seq1,
        composition=analyze_composition(seq1, record_id="seq1"),
        reverse_complement=reverse_complement(seq1),
        translation=translate_dna(seq1, record_id="seq1"),
    )

    seq2 = "GATTACA"
    r2 = RecordReport(
        record_id="seq2",
        description="seq2 Second test with\nnewline in header",
        sequence=seq2,
        composition=analyze_composition(seq2, record_id="seq2"),
        reverse_complement=reverse_complement(seq2),
        translation=translate_dna(seq2, record_id="seq2"),
    )

    return [r1, r2]


def test_format_terminal_report(sample_reports):
    """Test human-readable terminal report content."""
    text = format_terminal_report(sample_reports)
    assert "DNA SEQUENCE ANALYSIS TOOLKIT - REPORT" in text
    assert "[Record 1 of 2]: seq1" in text
    assert "GC-Content:       50.00%" in text
    assert "Rev Complement:   GCAT" in text
    assert "[Record 2 of 2]: seq2" in text
    assert "Length:           7 bp" in text


def test_export_json_structure(sample_reports):
    """Test JSON export produces valid JSON matching documented schema."""
    json_str = export_json(sample_reports, metadata={"test_run": True})
    data = json.loads(json_str)

    assert data["toolkit"] == "dna-toolkit"
    assert data["record_count"] == 2
    assert data["metadata"]["test_run"] is True
    assert len(data["records"]) == 2

    r0 = data["records"][0]
    assert r0["id"] == "seq1"
    assert r0["length"] == 4
    assert r0["composition"]["gc_percent"] == 50.0
    assert r0["reverse_complement"] == "GCAT"
    assert r0["translation"]["protein_sequence"] == "M"


def test_export_json_to_file(sample_reports, tmp_path: Path):
    """Test exporting JSON safely to a file."""
    out_file = tmp_path / "subdir" / "test_out.json"
    export_json(sample_reports, output_path=out_file)

    assert out_file.is_file()
    data = json.loads(out_file.read_text(encoding="utf-8"))
    assert data["record_count"] == 2


def test_export_csv_structure_and_rfc4180(sample_reports):
    """Test CSV export generates predictable columns and handles multiline descriptions cleanly."""
    csv_str = export_csv(sample_reports)
    reader = list(csv.DictReader(csv_str.splitlines()))

    assert len(reader) == 2
    expected_columns = {
        "record_id", "description", "length", "gc_count", "gc_percent",
        "count_A", "count_T", "count_G", "count_C", "count_ambiguous",
        "reverse_complement", "translation_frame", "translation_code", "protein_sequence"
    }
    assert set(reader[0].keys()) == expected_columns

    # Verify first row
    assert reader[0]["record_id"] == "seq1"
    assert reader[0]["length"] == "4"
    assert reader[0]["gc_percent"] == "50.0000"
    assert reader[0]["reverse_complement"] == "GCAT"

    # Verify second row (newlines in description replaced/normalized)
    assert "\n" not in reader[1]["description"]


def test_export_csv_to_file(sample_reports, tmp_path: Path):
    """Test exporting CSV safely to a file."""
    out_file = tmp_path / "output.csv"
    export_csv(sample_reports, output_path=out_file)

    assert out_file.is_file()
    lines = out_file.read_text(encoding="utf-8").splitlines()
    assert len(lines) == 3  # header + 2 records
