"""
Reporting and Export Module
===========================
Generates human-readable terminal reports, structured JSON exports, and
RFC-4180-compliant CSV exports for DNA analysis results.

Safe File Operations
--------------------
File writes utilize UTF-8 encoding with explicit error trapping for permission
and path errors to prevent data corruption.

CSV Escaping
------------
Per RFC 4180, fields containing delimiters, quotes, or newlines are escaped using
Python's standard `csv.writer`. Sequences and descriptions have line breaks stripped
to prevent corrupted tabular data.
"""

import csv
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import io
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from dna_toolkit.analysis import CompositionResult
from dna_toolkit.operations import TranslationResult


@dataclass
class RecordReport:
    """Consolidated analysis report for a single DNA sequence record."""
    record_id: str
    description: str
    sequence: str
    composition: Optional[CompositionResult] = None
    complement: Optional[str] = None
    reverse_complement: Optional[str] = None
    translation: Optional[TranslationResult] = None

    def to_dict(self) -> Dict[str, Any]:
        """Convert report to serializable dictionary."""
        data: Dict[str, Any] = {
            "id": self.record_id,
            "description": self.description,
            "length": len(self.sequence),
            "sequence": self.sequence,
        }
        if self.composition:
            data["composition"] = asdict(self.composition)
        if self.complement is not None:
            data["complement"] = self.complement
        if self.reverse_complement is not None:
            data["reverse_complement"] = self.reverse_complement
        if self.translation:
            data["translation"] = asdict(self.translation)
        return data


def format_terminal_report(reports: List[RecordReport]) -> str:
    """
    Format a list of RecordReport objects into a clean, human-readable terminal string.

    Args:
        reports: List of RecordReport objects.

    Returns:
        Formatted multi-line text report.
    """
    lines: List[str] = []
    divider = "=" * 70
    subdivider = "-" * 70

    lines.append(divider)
    lines.append(" DNA SEQUENCE ANALYSIS TOOLKIT - REPORT")
    lines.append(f" Total records analyzed: {len(reports)}")
    lines.append(divider)

    for i, report in enumerate(reports, start=1):
        lines.append(f"\n[Record {i} of {len(reports)}]: {report.record_id}")
        if report.description and report.description != report.record_id:
            lines.append(f" Description: {report.description}")
        lines.append(subdivider)

        # Sequence preview
        seq_display = report.sequence if len(report.sequence) <= 50 else f"{report.sequence[:47]}... (total {len(report.sequence)} bp)"
        lines.append(f" Sequence (5'->3'): {seq_display}")

        # Composition details
        if report.composition:
            comp = report.composition
            lines.append(f" Length:           {comp.length} bp")
            lines.append(f" GC-Content:       {comp.gc_percent:.2f}% (G+C: {comp.gc_count})")
            
            # Format counts cleanly
            base_counts_str = ", ".join(f"{b}: {c} ({comp.percentages.get(b, 0.0):.1f}%)" for b, c in comp.counts.items())
            lines.append(f" Nucleotides:      {base_counts_str}")
            if comp.ambiguous_count > 0:
                lines.append(f" Ambiguous Bases:  {comp.ambiguous_count} base(s)")
                lines.append(f" Unambiguous GC%:  {comp.unambiguous_gc_percent:.2f}%")

        # Operations
        if report.complement:
            comp_display = report.complement if len(report.complement) <= 50 else f"{report.complement[:47]}..."
            lines.append(f" Complement:       {comp_display} (3'->5')")
        if report.reverse_complement:
            rev_display = report.reverse_complement if len(report.reverse_complement) <= 50 else f"{report.reverse_complement[:47]}..."
            lines.append(f" Rev Complement:   {rev_display} (5'->3')")

        # Translation
        if report.translation:
            tr = report.translation
            prot_display = tr.protein_sequence if len(tr.protein_sequence) <= 50 else f"{tr.protein_sequence[:47]}... ({tr.codon_count} aa)"
            lines.append(f" Translation:      {prot_display}")
            lines.append(f"   [Frame: {tr.reading_frame:+d} | Code Table: {tr.genetic_code} | Codons: {tr.codon_count} | Stop Codon: {tr.has_stop_codon}]")
            if tr.partial_bases_trimmed > 0:
                lines.append(f"   * Partial trailing bases trimmed: {tr.partial_bases_trimmed} bp")

    lines.append("\n" + divider)
    lines.append(" END OF REPORT")
    lines.append(divider)

    return "\n".join(lines)


def export_json(
    reports: List[RecordReport],
    output_path: Optional[Union[str, Path]] = None,
    metadata: Optional[Dict[str, Any]] = None,
) -> str:
    """
    Export reports to a structured JSON string and optionally save to disk.

    Args:
        reports: List of RecordReport objects.
        output_path: Optional file path to write JSON file.
        metadata: Optional extra metadata dictionary.

    Returns:
        JSON string.
    """
    payload: Dict[str, Any] = {
        "toolkit": "dna-toolkit",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "record_count": len(reports),
        "metadata": metadata or {},
        "records": [r.to_dict() for r in reports],
    }

    json_str = json.dumps(payload, indent=2)

    if output_path:
        _safe_write_file(output_path, json_str)

    return json_str


def export_csv(
    reports: List[RecordReport],
    output_path: Optional[Union[str, Path]] = None,
) -> str:
    """
    Export reports to CSV format and optionally write to disk.

    Columns:
        record_id, description, length, gc_count, gc_percent, count_A, count_T,
        count_G, count_C, count_ambiguous, reverse_complement,
        translation_frame, translation_code, protein_sequence

    Args:
        reports: List of RecordReport objects.
        output_path: Optional file path to write CSV.

    Returns:
        CSV string content.
    """
    buffer = io.StringIO()
    fieldnames = [
        "record_id",
        "description",
        "length",
        "gc_count",
        "gc_percent",
        "count_A",
        "count_T",
        "count_G",
        "count_C",
        "count_ambiguous",
        "reverse_complement",
        "translation_frame",
        "translation_code",
        "protein_sequence",
    ]

    writer = csv.DictWriter(buffer, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()

    for r in reports:
        row: Dict[str, Any] = {
            "record_id": r.record_id,
            "description": " ".join(r.description.split()),  # normalize newlines
            "length": len(r.sequence),
            "gc_count": r.composition.gc_count if r.composition else "",
            "gc_percent": f"{r.composition.gc_percent:.4f}" if r.composition else "",
            "count_A": r.composition.counts.get("A", "") if r.composition else "",
            "count_T": r.composition.counts.get("T", "") if r.composition else "",
            "count_G": r.composition.counts.get("G", "") if r.composition else "",
            "count_C": r.composition.counts.get("C", "") if r.composition else "",
            "count_ambiguous": r.composition.ambiguous_count if r.composition else "",
            "reverse_complement": r.reverse_complement or "",
            "translation_frame": r.translation.reading_frame if r.translation else "",
            "translation_code": r.translation.genetic_code if r.translation else "",
            "protein_sequence": r.translation.protein_sequence if r.translation else "",
        }
        writer.writerow(row)

    csv_str = buffer.getvalue()

    if output_path:
        _safe_write_file(output_path, csv_str)

    return csv_str


def _safe_write_file(file_path: Union[str, Path], content: str) -> None:
    """
    Safely write string content to disk with directory creation and explicit error trapping.

    Args:
        file_path: Destination path.
        content: Text content to write.

    Raises:
        OSError: If destination cannot be written.
    """
    target = Path(file_path)
    try:
        target.parent.mkdir(parents=True, exist_ok=True)
        with open(target, "w", encoding="utf-8") as f:
            f.write(content)
    except PermissionError as err:
        raise PermissionError(f"Permission denied writing to output file: {target}") from err
    except OSError as err:
        raise OSError(f"Failed to write output file at '{target}': {err}") from err
