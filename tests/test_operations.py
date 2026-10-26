"""
Tests for Sequence Operations and Translation
=============================================
Validates complement vs reverse complement, 6-frame translation, stop codons,
incomplete codons, and alternative genetic codes.
"""

import pytest

from dna_toolkit.operations import (
    complement,
    reverse_complement,
    translate_dna,
)


def test_complement_vs_reverse_complement_distinction():
    """
    Verify known case:
      AGTC complement is TCAG (3' to 5').
      AGTC reverse complement is GACT (5' to 3').
    """
    seq = "AGTC"
    comp = complement(seq)
    rev_comp = reverse_complement(seq)

    assert comp == "TCAG"
    assert rev_comp == "GACT"
    assert comp != rev_comp


def test_known_translation_atggaa():
    """Verify known case: ATGGAA translates to ME in standard genetic code."""
    result = translate_dna("ATGGAA")
    assert result.protein_sequence == "ME"
    assert result.codon_count == 2
    assert result.has_stop_codon is False
    assert result.partial_bases_trimmed == 0


def test_six_reading_frames():
    """
    Test forward and reverse reading frames.
    Sequence: ATG GCT AGC TAA
    Length: 12 bp
    Frame +1: ATG GCT AGC TAA -> M A S *
    Frame +2: TGG CTA GCT AA  -> W L A (trimmed 2 bases: 'AA')
    Frame +3: GGC TAG CTA A   -> G * L (trimmed 1 base: 'A')
    """
    seq = "ATGGCTAGCTAA"

    # Frame +1
    res_f1 = translate_dna(seq, frame=1)
    assert res_f1.protein_sequence == "MAS*"
    assert res_f1.has_stop_codon is True
    assert res_f1.partial_bases_trimmed == 0

    # Frame +2
    res_f2 = translate_dna(seq, frame=2)
    assert res_f2.protein_sequence == "WLA"
    assert res_f2.partial_bases_trimmed == 2

    # Frame +3
    res_f3 = translate_dna(seq, frame=3)
    assert res_f3.protein_sequence == "G*L"
    assert res_f3.partial_bases_trimmed == 1

    # Reverse frame -1 (Reverse complement of ATGGCTAGCTAA is TTAGCTAGCCAT)
    # TTA GCT AGC CAT -> L A S H
    res_r1 = translate_dna(seq, frame=-1)
    assert res_r1.protein_sequence == "LASH"


def test_translation_to_stop():
    """Test halting translation at the first stop codon when to_stop=True."""
    seq = "ATGGCTTAAGCC"  # ATG (M), GCT (A), TAA (*), GCC (A)
    # Default (to_stop=False): translates through stop codon
    full_res = translate_dna(seq, to_stop=False)
    assert full_res.protein_sequence == "MA*A"

    # to_stop=True: stops at TAA without including '*'
    stop_res = translate_dna(seq, to_stop=True)
    assert stop_res.protein_sequence == "MA"
    assert stop_res.codon_count == 2


def test_incomplete_codon_handling():
    """Verify that trailing 1 or 2 bases are cleanly trimmed with recorded metadata."""
    # Length 7: 2 codons + 1 base
    result = translate_dna("ATGAAAT")
    assert result.protein_sequence == "MK"
    assert result.codon_count == 2
    assert result.partial_bases_trimmed == 1

    # Length 8: 2 codons + 2 bases
    result2 = translate_dna("ATGAAATG")
    assert result2.protein_sequence == "MK"
    assert result2.partial_bases_trimmed == 2


def test_alternative_genetic_code():
    """
    Verify translation using NCBI alternative genetic code table.
    Table 2: Vertebrate Mitochondrial Code
    In Standard Code (Table 1):
      AGA is Arginine ('R')
      UGA/TGA is Stop ('*')
    In Vertebrate Mitochondrial Code (Table 2):
      AGA is Stop ('*')
      TGA is Tryptophan ('W')
    """
    seq = "TGATGA"  # Two TGA codons
    std_res = translate_dna(seq, table=1)
    assert std_res.protein_sequence == "**"

    mito_res = translate_dna(seq, table=2)
    assert mito_res.protein_sequence == "WW"


def test_invalid_reading_frame_raises_value_error():
    """Verify that invalid frames (e.g. 0, 4, -4) raise ValueError."""
    with pytest.raises(ValueError) as exc_info:
        translate_dna("ATGC", frame=4)
    assert "Invalid reading frame" in str(exc_info.value)

    with pytest.raises(ValueError):
        translate_dna("ATGC", frame=0)
