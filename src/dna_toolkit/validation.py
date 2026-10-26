"""
DNA Sequence Validation Module
==============================
Provides explicit, rule-based validation policies for DNA sequence data.

Why Sequence Validation is Critical in Bioinformatics
------------------------------------------------------
Real-world genomic data is frequently noisy. Raw files can contain:
  1. Non-nucleotide characters (e.g., numbers from numbered FASTA, asterisks from
     translated proteins mistakenly passed as DNA, punctuation, or spaces).
  2. Degenerate/Ambiguous IUPAC codes (N for any base, R for purine, Y for pyrimidine).
     These arise from sequencer uncertainty (low Phred scores) or natural genetic
     heterozygosity.
  3. RNA symbols (Uracil 'U' instead of Thymine 'T').
  4. Empty sequence records where headers exist without following nucleotide lines.

If unvalidated data is passed to downstream pipelines, calculations like GC-content
produce misleading denominators, and translation tables will fail or translate invalid
codons into unintended amino acids.

Validation Policies
-------------------
- STRICT: Only standard canonical deoxyribonucleotides (A, C, G, T) are accepted.
- AMBIGUOUS: Standard canonical nucleotides plus IUPAC degenerate base codes are accepted:
    - R: A or G (puRine)
    - Y: C or T (pYrimidine)
    - S: G or C (Strong - 3 H-bonds)
    - W: A or T (Weak - 2 H-bonds)
    - K: G or T (Keto)
    - M: A or C (aMino)
    - B: C or G or T (not A)
    - D: A or G or T (not C)
    - H: A or C or T (not G)
    - V: A or C or G (not T)
    - N: Any base (A, C, G, or T)

Duplicate Identifier Policy
---------------------------
In bioinformatics, FASTA headers define sequence namespaces. Duplicate identifiers
compromise data provenance and downstream indexing. By default, `validate_records`
enforces strict identifier uniqueness and raises `DuplicateRecordError`. If
`allow_duplicates=True` is explicitly passed, duplicates are permitted and documented.
"""

from enum import Enum
from typing import Dict, List, Set, Union

from dna_toolkit.parser import DnaRecord

# Canonical DNA bases
CANONICAL_BASES: Set[str] = {"A", "C", "G", "T"}

# IUPAC degenerate nucleotide codes for DNA (includes canonical)
IUPAC_AMBIGUOUS_BASES: Set[str] = {
    "A", "C", "G", "T",  # Canonical
    "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N",  # Degenerate
}


class ValidationMode(str, Enum):
    """Supported validation policies."""
    STRICT = "strict"
    AMBIGUOUS = "ambiguous"


class ValidationError(Exception):
    """Base exception for all sequence validation failures."""
    pass


class EmptySequenceError(ValidationError):
    """Raised when a record contains an empty sequence (0 nucleotides)."""
    pass


class InvalidBaseError(ValidationError):
    """
    Raised when invalid sequence characters are encountered.

    Attributes:
        invalid_chars: Set of invalid characters discovered.
        char_positions: Mapping of invalid character to 1-based index positions.
    """
    def __init__(self, message: str, invalid_chars: Set[str], char_positions: Dict[str, List[int]]):
        super().__init__(message)
        self.invalid_chars = invalid_chars
        self.char_positions = char_positions


class DuplicateRecordError(ValidationError):
    """Raised when duplicate sequence identifiers are detected in a FASTA dataset."""
    pass


def validate_sequence(
    sequence: str,
    mode: Union[ValidationMode, str] = ValidationMode.STRICT,
    record_id: str = "unnamed_record",
) -> str:
    """
    Validate a single DNA sequence string against the specified policy.

    Case normalization:
        Lower-case characters (e.g., 'atgc') are automatically converted to uppercase ('ATGC').
        In genomics, lower-case is often used to denote repeat-masked regions (e.g., RepeatMasker),
        which are biologically valid DNA.

    Args:
        sequence: The raw DNA sequence string.
        mode: ValidationMode.STRICT or ValidationMode.AMBIGUOUS.
        record_id: Optional record identifier for actionable error reporting.

    Returns:
        The normalized, uppercase DNA sequence string.

    Raises:
        EmptySequenceError: If sequence length is 0.
        InvalidBaseError: If symbols outside the allowed set are found.
    """
    mode_enum = ValidationMode(mode) if isinstance(mode, str) else mode
    cleaned = "".join(sequence.split()).upper()

    if not cleaned:
        raise EmptySequenceError(
            f"Record '{record_id}' contains an empty sequence (0 nucleotides)."
        )

    allowed_set = CANONICAL_BASES if mode_enum == ValidationMode.STRICT else IUPAC_AMBIGUOUS_BASES

    invalid_positions: Dict[str, List[int]] = {}
    for idx, char in enumerate(cleaned, start=1):
        if char not in allowed_set:
            invalid_positions.setdefault(char, []).append(idx)

    if invalid_positions:
        bad_chars = sorted(invalid_positions.keys())
        examples = []
        for c in bad_chars[:5]:
            positions = invalid_positions[c][:3]
            pos_str = ", ".join(str(p) for p in positions)
            if len(invalid_positions[c]) > 3:
                pos_str += ", ..."
            examples.append(f"'{c}' at position(s) {pos_str}")
        
        detail = "; ".join(examples)
        policy_desc = "canonical DNA (A, T, G, C)" if mode_enum == ValidationMode.STRICT else "IUPAC DNA (A, T, G, C + ambiguous codes)"
        
        raise InvalidBaseError(
            f"Record '{record_id}' failed {mode_enum.value} validation. "
            f"Found invalid characters for {policy_desc}: {detail}.",
            invalid_chars=set(invalid_positions.keys()),
            char_positions=invalid_positions,
        )

    return cleaned


def validate_records(
    records: List[DnaRecord],
    mode: Union[ValidationMode, str] = ValidationMode.STRICT,
    allow_duplicates: bool = False,
) -> List[DnaRecord]:
    """
    Validate a collection of DnaRecord objects.

    Args:
        records: List of parsed DnaRecord objects.
        mode: Validation policy (STRICT or AMBIGUOUS).
        allow_duplicates: If False (default), duplicate record IDs raise DuplicateRecordError.

    Returns:
        A list of new DnaRecord objects containing normalized, uppercase sequences.

    Raises:
        DuplicateRecordError: If duplicate sequence IDs are found and allow_duplicates is False.
        EmptySequenceError: If any sequence is empty.
        InvalidBaseError: If any sequence contains unauthorized characters.
    """
    seen_ids: Set[str] = set()
    validated_records: List[DnaRecord] = []

    for record in records:
        if not allow_duplicates:
            if record.id in seen_ids:
                raise DuplicateRecordError(
                    f"Duplicate sequence identifier '{record.id}' detected in dataset. "
                    "Identifiers must be unique unless --allow-duplicates is enabled."
                )
            seen_ids.add(record.id)

        valid_seq = validate_sequence(record.sequence, mode=mode, record_id=record.id)
        validated_records.append(
            DnaRecord(
                id=record.id,
                description=record.description,
                sequence=valid_seq,
            )
        )

    return validated_records
