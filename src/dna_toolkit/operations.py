"""
DNA Sequence Operations Module
==============================
Provides biological sequence transformations:
  - Strand complementation and reverse-complementation.
  - DNA-to-protein 6-frame translation using Biopython.

Biological Concepts and Decisions
---------------------------------

1. Complement vs. Reverse Complement:
   DNA double helix consists of two antiparallel strands:
     Forward strand:  5'- A G T C - 3'
     Direct complement: 3'- T C A G - 5'
   In biology and bioinformatics, nucleic acids are universally written 5' to 3'.
   Therefore, to read the opposite strand in its proper 5' to 3' direction, we reverse
   the complement:
     Reverse complement: 5'- G A C T - 3'
   It is vital never to confuse direct complement with reverse complement!

2. Why Translation Requires Reading Frames (6-frame translation):
   An arbitrary DNA sequence read from a FASTA file is merely a fragment of a double-stranded
   chromosome. We cannot assume:
     a) That the sequence starts with an AUG/ATG start codon.
     b) Which reading frame the organism's ribosome uses.
     c) Whether the protein-coding gene resides on the forward or reverse strand.

   Double-stranded DNA therefore has six potential open reading frames (ORFs):
     - Forward Strand:
         Frame +1: Begins at index 0 (bases 0, 1, 2 form codon 1)
         Frame +2: Begins at index 1 (bases 1, 2, 3 form codon 1)
         Frame +3: Begins at index 2 (bases 2, 3, 4 form codon 1)
     - Reverse Complement Strand:
         Frame -1: Begins at index 0 of reverse complement
         Frame -2: Begins at index 1 of reverse complement
         Frame -3: Begins at index 2 of reverse complement

3. Incomplete and Ambiguous Codons:
   A complete codon requires exactly 3 nucleotides. If a reading frame's length is not
   a multiple of 3, the trailing 1 or 2 bases cannot form an amino acid.
   Rather than letting translation fail or emit unhandled warnings, we explicitly trim
   and record partial trailing bases (`partial_bases_trimmed`).

4. Stop Codons:
   Standard stop codons (TAA, TAG, TGA) are represented by an asterisk ('*').
   When `to_stop=True`, translation halts at the first encountered stop codon.
"""

from dataclasses import dataclass
from typing import Union

from Bio.Seq import Seq

from dna_toolkit.validation import (
    EmptySequenceError,
    ValidationMode,
    validate_sequence,
)


@dataclass(frozen=True)
class TranslationResult:
    """
    Structured outcome of DNA-to-protein translation.

    Attributes:
        record_id: Identifier of the analyzed record.
        reading_frame: Reading frame used (+1, +2, +3, -1, -2, -3).
        genetic_code: NCBI genetic code table ID or name (default 1 / 'Standard').
        protein_sequence: Translated amino acid single-letter sequence.
        codon_count: Number of codons successfully translated.
        has_stop_codon: True if one or more stop codons ('*') are present.
        partial_bases_trimmed: Number of incomplete trailing nucleotides trimmed (0, 1, or 2).
        to_stop: Whether translation terminated at the first stop codon.
    """
    record_id: str
    reading_frame: int
    genetic_code: Union[str, int]
    protein_sequence: str
    codon_count: int
    has_stop_codon: bool
    partial_bases_trimmed: int
    to_stop: bool


def complement(sequence: str, mode: Union[ValidationMode, str] = ValidationMode.STRICT) -> str:
    """
    Generate the direct Watson-Crick complement of a DNA sequence (3' to 5').

    Example:
        'AGTC' -> 'TCAG'

    Args:
        sequence: DNA sequence string.
        mode: Validation policy.

    Returns:
        Complemented sequence string.
    """
    valid_seq = validate_sequence(sequence, mode=mode)
    bio_seq = Seq(valid_seq)
    return str(bio_seq.complement())


def reverse_complement(sequence: str, mode: Union[ValidationMode, str] = ValidationMode.STRICT) -> str:
    """
    Generate the reverse complement of a DNA sequence (5' to 3' of opposite strand).

    Example:
        'AGTC' -> 'GACT'

    Args:
        sequence: DNA sequence string.
        mode: Validation policy.

    Returns:
        Reverse complement sequence string.
    """
    valid_seq = validate_sequence(sequence, mode=mode)
    bio_seq = Seq(valid_seq)
    return str(bio_seq.reverse_complement())


def translate_dna(
    sequence: str,
    frame: int = 1,
    table: Union[str, int] = 1,
    to_stop: bool = False,
    mode: Union[ValidationMode, str] = ValidationMode.STRICT,
    record_id: str = "unnamed_record",
) -> TranslationResult:
    """
    Translate a DNA sequence to a protein sequence using Biopython.

    Args:
        sequence: DNA sequence string.
        frame: Reading frame (+1, +2, +3 for forward; -1, -2, -3 for reverse strand).
        table: NCBI genetic code table (default 1: 'Standard').
        to_stop: If True, stop translation at the first stop codon (excluding '*').
        mode: Validation policy (STRICT or AMBIGUOUS).
        record_id: Identifier of the record being translated.

    Returns:
        TranslationResult object containing the translated protein and metadata.

    Raises:
        ValueError: If frame is not in {1, 2, 3, -1, -2, -3}.
        EmptySequenceError: If sequence is empty.
        InvalidBaseError: If sequence fails validation.
    """
    valid_frames = {1, 2, 3, -1, -2, -3}
    if frame not in valid_frames:
        raise ValueError(
            f"Invalid reading frame: {frame}. Must be one of {sorted(valid_frames)}."
        )

    valid_seq = validate_sequence(sequence, mode=mode, record_id=record_id)
    if not valid_seq:
        raise EmptySequenceError(f"Record '{record_id}' has empty sequence.")

    # Determine working sequence based on strand
    if frame > 0:
        working_seq = valid_seq
        offset = frame - 1
    else:
        working_seq = str(Seq(valid_seq).reverse_complement())
        offset = abs(frame) - 1

    # Extract reading frame subsequence
    subseq = working_seq[offset:]

    # Handle partial trailing codons explicitly
    remainder = len(subseq) % 3
    if remainder != 0:
        complete_subseq = subseq[:-remainder]
    else:
        complete_subseq = subseq

    if not complete_subseq:
        # If sequence was too short for even 1 complete codon in this frame
        protein = ""
        codon_count = 0
        has_stop = False
    else:
        bio_seq = Seq(complete_subseq)
        protein = str(bio_seq.translate(table=table, to_stop=to_stop))
        codon_count = len(protein)
        has_stop = "*" in protein

    return TranslationResult(
        record_id=record_id,
        reading_frame=frame,
        genetic_code=table,
        protein_sequence=protein,
        codon_count=codon_count,
        has_stop_codon=has_stop,
        partial_bases_trimmed=remainder,
        to_stop=to_stop,
    )
