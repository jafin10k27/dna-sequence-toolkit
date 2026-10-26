"""
FASTA Parser Module
===================
Handles reading and parsing biological sequence records from FASTA-formatted files
using Biopython's SeqIO engine, complemented with explicit sanity checks.

Why Biopython?
--------------
FASTA is a standard file format in bioinformatics dating back to 1985 (David J. Lipman
and William R. Pearson). While conceptually simple (lines starting with '>' contain
metadata, followed by lines of sequence characters), real-world FASTA files have
subtleties:
  1. Multiline sequences: A single 3-billion-base human chromosome sequence is wrapped
     across millions of lines (typically 60-80 chars per line).
  2. Whitespace and line endings: Windows (CRLF), Unix (LF), and legacy Mac (CR) line
     endings must be handled portably.
  3. Memory efficiency: Biopython uses generators (`SeqIO.parse`) to stream records
     without loading an entire multi-gigabyte FASTA file into RAM.

In this module, we wrap `Bio.SeqIO.parse` with strict checks for missing files, empty
files, and malformed content so the user receives actionable error messages.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import List, Union

from Bio import SeqIO


class FastaError(Exception):
    """Base exception for all FASTA parsing errors."""
    pass


class FastaFileNotFoundError(FastaError, FileNotFoundError):
    """Raised when the specified FASTA file does not exist."""
    pass


class EmptyFastaError(FastaError):
    """Raised when the FASTA file exists but contains no data (0 bytes or whitespace)."""
    pass


class MalformedFastaError(FastaError):
    """Raised when a file cannot be parsed as valid FASTA (e.g., missing '>' headers)."""
    pass


@dataclass(frozen=True)
class DnaRecord:
    """
    Immutable representation of an individual biological sequence record.

    Attributes:
        id: Unique or primary identifier for the sequence (e.g., 'seq1', 'NM_000518.5').
        description: Full header line information, including annotations and titles.
        sequence: Raw sequence string as read from the FASTA file (not yet validated).
    """
    id: str
    description: str
    sequence: str

    def __post_init__(self):
        # Ensure sequence characters are preserved exactly as read (whitespace stripped)
        object.__setattr__(self, "sequence", self.sequence.strip())
        object.__setattr__(self, "id", self.id.strip())
        object.__setattr__(self, "description", self.description.strip())


def parse_fasta(file_path: Union[str, Path]) -> List[DnaRecord]:
    """
    Parse a FASTA file and return a list of DnaRecord objects.

    Args:
        file_path: Path to the FASTA file (string or Path object).

    Returns:
        List of DnaRecord instances containing id, description, and sequence.

    Raises:
        FastaFileNotFoundError: If the file does not exist on disk.
        EmptyFastaError: If the file is 0 bytes or contains only whitespace.
        MalformedFastaError: If the file does not contain valid FASTA records (e.g., no '>' headers).
    """
    path = Path(file_path)

    if not path.is_file():
        raise FastaFileNotFoundError(f"FASTA file not found at path: {path.resolve()}")

    # Check file size before opening
    if path.stat().st_size == 0:
        raise EmptyFastaError(f"FASTA file is empty (0 bytes): {path.name}")

    # Inspect file content for empty/whitespace-only files or lack of FASTA header
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        first_chunk = handle.read(4096)
        if not first_chunk.strip():
            raise EmptyFastaError(f"FASTA file contains only whitespace: {path.name}")
        
        # In FASTA format, valid files must begin with '>' on the first non-empty line
        first_non_empty = next((line.strip() for line in first_chunk.splitlines() if line.strip()), "")
        if not first_non_empty.startswith(">"):
            raise MalformedFastaError(
                f"File '{path.name}' does not appear to be a FASTA file. "
                f"Expected header starting with '>', found: '{first_non_empty[:30]}...'"
            )

    records: List[DnaRecord] = []
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            for seq_record in SeqIO.parse(handle, "fasta"):
                raw_seq = str(seq_record.seq)
                records.append(
                    DnaRecord(
                        id=str(seq_record.id),
                        description=str(seq_record.description),
                        sequence=raw_seq,
                    )
                )
    except Exception as exc:
        raise MalformedFastaError(f"Failed to parse FASTA file '{path.name}': {exc}") from exc

    if not records:
        raise MalformedFastaError(f"No valid FASTA records found in '{path.name}'.")

    return records
