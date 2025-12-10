#!/usr/bin/env python3
"""Pytest test suite for protein_metrics.py"""

import pytest
from protein_metrics import (
    _clean_sequence,
    length,
    amino_acid_composition,
    molecular_weight,
    get_summary_statistics,
    get_amino_acid_percentage,
    ProteinMetricsError,
    AMINO_ACID_WEIGHTS
)


class TestCleanSequence:
    """Tests for the _clean_sequence helper function."""
    
    def test_uppercase_conversion(self):
        assert _clean_sequence("mklt") == "MKLT"
    
    def test_whitespace_removal(self):
        assert _clean_sequence("MKLT VVGA\nGGKS") == "MKLTVVGAGGKS"
    
    def test_mixed_case_and_whitespace(self):
        assert _clean_sequence("  mklt vvga  ") == "MKLTVVGA"


class TestLength:
    """Tests for length function."""
    
    def test_simple_sequence(self):
        assert length("MKLVVVGA") == 8
    
    def test_with_whitespace(self):
        assert length("MKLT VVGA\nGGKS") == 12
    
    def test_empty_sequence_raises(self):
        with pytest.raises(ProteinMetricsError, match="empty or contains no valid characters"):
            length("")
    
    def test_only_whitespace_raises(self):
        with pytest.raises(ProteinMetricsError):
            length("   \n\t  ")


class TestAminoAcidComposition:
    """Tests for amino_acid_composition function."""
    
    def test_simple_composition(self):
        comp = amino_acid_composition("MKLVVVGA")
        assert comp['M'] == 1
        assert comp['K'] == 1
        assert comp['L'] == 1
        assert comp['V'] == 3
        assert comp['G'] == 1
        assert comp['A'] == 1
        assert len(comp) == 6  # 6 unique amino acids
        assert sum(comp.values()) == 8  # Total residues
    
    def test_single_amino_acid(self):
        comp = amino_acid_composition("AAA")
        assert comp == {'A': 3}
    
    def test_all_20_amino_acids(self):
        comp = amino_acid_composition("ACDEFGHIKLMNPQRSTVWY")
        assert len(comp) == 20
        for aa in comp:
            assert comp[aa] == 1
    
    def test_with_whitespace(self):
        comp1 = amino_acid_composition("MKLVVVGA")
        comp2 = amino_acid_composition("MKL VVV\nGA")
        assert comp1 == comp2
    
    def test_case_insensitive(self):
        comp1 = amino_acid_composition("MKLVVVGA")
        comp2 = amino_acid_composition("mklvvvga")
        assert comp1 == comp2
    
    def test_empty_sequence_raises(self):
        with pytest.raises(ProteinMetricsError):
            amino_acid_composition("")
    
    def test_no_valid_amino_acids_raises(self):
        # "123" has no alphabetic characters after cleaning
        with pytest.raises(ProteinMetricsError, match="empty or contains no valid characters"):
            amino_acid_composition("123")
    
    def test_only_non_standard_raises(self):
        # X, B, Z, J are not in the standard 20
        with pytest.raises(ProteinMetricsError):
            amino_acid_composition("XXXX")


class TestMolecularWeight:
    """Tests for molecular_weight function."""
    
    def test_three_alanines(self):
        # 3 * 89.09 + 18.015 = 285.285
        mw = molecular_weight("AAA")
        expected = 3 * AMINO_ACID_WEIGHTS['A'] + 18.015
        assert mw == pytest.approx(expected, rel=0.001)
    
    def test_three_glycines(self):
        # 3 * 75.07 + 18.015 = 243.225
        mw = molecular_weight("GGG")
        expected = 3 * AMINO_ACID_WEIGHTS['G'] + 18.015
        assert mw == pytest.approx(expected, rel=0.001)
    
    def test_mixed_sequence(self):
        # Calculate expected weight manually
        seq = "MKLT"
        expected = sum(AMINO_ACID_WEIGHTS[aa] for aa in seq) + 18.015
        mw = molecular_weight(seq)
        assert mw == pytest.approx(expected, rel=0.001)
    
    def test_with_whitespace(self):
        mw1 = molecular_weight("MKLVVVGA")
        mw2 = molecular_weight("MKL VVV\nGA")
        assert mw1 == pytest.approx(mw2)
    
    def test_empty_sequence_raises(self):
        with pytest.raises(ProteinMetricsError):
            molecular_weight("")
    
    def test_no_valid_amino_acids_raises(self):
        with pytest.raises(ProteinMetricsError):
            molecular_weight("123")


class TestGetSummaryStatistics:
    """Tests for get_summary_statistics function."""
    
    def test_comprehensive_stats(self):
        seq = "MKLVVVGA"
        stats = get_summary_statistics(seq)
        
        # Verify all keys are present
        assert 'length' in stats
        assert 'molecular_weight' in stats
        assert 'composition' in stats
        assert 'unique_amino_acids' in stats
        
        # Verify values match individual function calls
        assert stats['length'] == length(seq)
        assert stats['molecular_weight'] == pytest.approx(molecular_weight(seq))
        assert stats['composition'] == amino_acid_composition(seq)
        assert stats['unique_amino_acids'] == len(amino_acid_composition(seq))
    
    def test_all_20_amino_acids(self):
        seq = "ACDEFGHIKLMNPQRSTVWY"
        stats = get_summary_statistics(seq)
        
        assert stats['length'] == 20
        assert stats['unique_amino_acids'] == 20
        assert len(stats['composition']) == 20
    
    def test_simple_peptide(self):
        seq = "MKLVVVGAGGKS"
        stats = get_summary_statistics(seq)
        
        assert stats['length'] == 12
        assert stats['molecular_weight'] > 0
        assert sum(stats['composition'].values()) == 12


class TestGetAminoAcidPercentage:
    """Tests for get_amino_acid_percentage function."""
    
    def test_percentages_sum_to_100(self):
        seq = "AAAGGG"
        pct = get_amino_acid_percentage(seq)
        assert sum(pct.values()) == pytest.approx(100.0)
    
    def test_equal_distribution(self):
        seq = "AAAGGG"
        pct = get_amino_acid_percentage(seq)
        assert pct['A'] == pytest.approx(50.0)
        assert pct['G'] == pytest.approx(50.0)
    
    def test_single_amino_acid(self):
        seq = "AAA"
        pct = get_amino_acid_percentage(seq)
        assert pct == {'A': pytest.approx(100.0)}
    
    def test_all_20_amino_acids(self):
        seq = "ACDEFGHIKLMNPQRSTVWY"
        pct = get_amino_acid_percentage(seq)
        # Each amino acid should be 5%
        for aa_pct in pct.values():
            assert aa_pct == pytest.approx(5.0)
    
    def test_empty_sequence_raises(self):
        with pytest.raises(ProteinMetricsError):
            get_amino_acid_percentage("")


class TestEdgeCases:
    """Tests for various edge cases."""
    
    def test_very_long_sequence(self):
        # Test with a long sequence
        seq = "A" * 1000
        assert length(seq) == 1000
        assert amino_acid_composition(seq) == {'A': 1000}
    
    def test_single_residue(self):
        seq = "M"
        comp = amino_acid_composition(seq)
        assert comp == {'M': 1}
        assert length(seq) == 1
        mw = molecular_weight(seq)
        assert mw == pytest.approx(AMINO_ACID_WEIGHTS['M'] + 18.015)
    
    def test_glycine_rich(self):
        # Glycine is the lightest amino acid
        seq = "GGGGGGGGGG"
        mw = molecular_weight(seq)
        expected = 10 * AMINO_ACID_WEIGHTS['G'] + 18.015
        assert mw == pytest.approx(expected)
    
    def test_tryptophan_rich(self):
        # Tryptophan is the heaviest amino acid
        seq = "WWWWW"
        mw = molecular_weight(seq)
        expected = 5 * AMINO_ACID_WEIGHTS['W'] + 18.015
        assert mw == pytest.approx(expected)
