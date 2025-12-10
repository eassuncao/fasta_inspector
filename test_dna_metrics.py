#!/usr/bin/env python3
"""Pytest test suite for dna_metrics.py"""

import pytest
from dna_metrics import (
    _clean_sequence,
    length,
    gc_content,
    n_count,
    at_gc_ratio,
    get_base_counts,
    get_summary_statistics
)


class TestCleanSequence:
    """Tests for the _clean_sequence helper function."""
    
    def test_uppercase_conversion(self):
        assert _clean_sequence("atcg") == "ATCG"
    
    def test_whitespace_removal(self):
        assert _clean_sequence("ATCG ATCG\nATCG") == "ATCGATCGATCG"
    
    def test_mixed_case_and_whitespace(self):
        assert _clean_sequence("  atcg ATCG\n  ") == "ATCGATCG"


class TestLength:
    """Tests for length function."""
    
    def test_simple_sequence(self):
        assert length("ATCGATCG") == 8
    
    def test_with_whitespace(self):
        assert length("ATCG ATCG\nATCG") == 12
    
    def test_empty_sequence(self):
        assert length("") == 0
    
    def test_only_whitespace(self):
        assert length("   \n\t  ") == 0


class TestGCContent:
    """Tests for gc_content function."""
    
    def test_balanced_dna(self):
        assert gc_content("ATCGATCG") == pytest.approx(50.0)
    
    def test_gc_rich(self):
        assert gc_content("GGGGCCCC") == pytest.approx(100.0)
    
    def test_at_rich(self):
        assert gc_content("AAAAATTTTT") == pytest.approx(0.0)
    
    def test_with_ambiguity(self):
        # ATCGNNNATCG: 11 total bases, 4 are G or C
        assert gc_content("ATCGNNNATCG") == pytest.approx(36.36, rel=0.01)
    
    def test_empty_sequence(self):
        assert gc_content("") == 0.0
    
    def test_only_n(self):
        assert gc_content("NNNNNN") == 0.0
    
    def test_with_whitespace(self):
        assert gc_content("ATCG ATCG\nATCG") == pytest.approx(50.0)
    
    def test_rna_uracil(self):
        # AUCGAUCG: U is not counted in GC, but in total length
        assert gc_content("AUCGAUCG") == pytest.approx(50.0)


class TestNCount:
    """Tests for n_count function."""
    
    def test_no_ambiguity(self):
        assert n_count("ATCGATCG") == 0
    
    def test_with_ns(self):
        assert n_count("ATCGNNNATCG") == 3
    
    def test_only_ns(self):
        assert n_count("NNNNNN") == 6
    
    def test_empty_sequence(self):
        assert n_count("") == 0
    
    def test_case_insensitive(self):
        assert n_count("ATCGnnnATCG") == 3


class TestATGCRatio:
    """Tests for at_gc_ratio function."""
    
    def test_balanced(self):
        # ATCG: 2 AT, 2 GC -> ratio 1.0
        assert at_gc_ratio("ATCG") == pytest.approx(1.0)
    
    def test_at_rich_ratio(self):
        # AAATCG: 4 AT, 2 GC -> ratio 2.0
        assert at_gc_ratio("AAATCG") == pytest.approx(2.0)
    
    def test_only_at_no_gc(self):
        # Only AT bases -> infinity
        assert at_gc_ratio("AAAAATTTTT") == float('inf')
    
    def test_only_gc_no_at(self):
        # Only GC bases -> 0.0
        assert at_gc_ratio("GGGGCCCC") == 0.0
    
    def test_empty_sequence(self):
        assert at_gc_ratio("") == 0.0
    
    def test_only_n(self):
        # Only N bases -> 0 AT, 0 GC -> 0.0
        assert at_gc_ratio("NNNNNN") == 0.0
    
    def test_with_ambiguity(self):
        # ATCGNNNATCG: 4 AT, 4 GC -> ratio 1.0
        assert at_gc_ratio("ATCGNNNATCG") == pytest.approx(1.0)
    
    def test_rna_uracil_counted_as_t(self):
        # AUCG: A and U (counted as T) = 2 AT, C and G = 2 GC -> ratio 1.0
        assert at_gc_ratio("AUCG") == pytest.approx(1.0)
    
    def test_mixed_rna_dna(self):
        # ATCGUUUATCG: A=2, T=2, U=3 (as T), C=2, G=2
        # AT total = 2+2+3 = 7, GC total = 4 -> ratio 1.75
        assert at_gc_ratio("ATCGUUUATCG") == pytest.approx(1.75)


class TestGetBaseCounts:
    """Tests for get_base_counts function."""
    
    def test_simple_dna(self):
        counts = get_base_counts("ATCGATCG")
        assert counts['A'] == 2
        assert counts['T'] == 2
        assert counts['C'] == 2
        assert counts['G'] == 2
        assert counts['N'] == 0
        assert counts['other'] == 0
        assert counts['total'] == 8
    
    def test_with_ns(self):
        counts = get_base_counts("ATCGNNNATCG")
        assert counts['A'] == 2
        assert counts['T'] == 2
        assert counts['C'] == 2
        assert counts['G'] == 2
        assert counts['N'] == 3
        assert counts['other'] == 0
        assert counts['total'] == 11
    
    def test_rna_uracil_counted_as_t(self):
        counts = get_base_counts("AUCGAUCG")
        assert counts['A'] == 2
        assert counts['T'] == 2  # U counted as T
        assert counts['C'] == 2
        assert counts['G'] == 2
        assert counts['N'] == 0
        assert counts['total'] == 8
    
    def test_with_iupac_ambiguity_codes(self):
        # ATCGRYSWKM: R,Y,S,W,K,M are IUPAC codes -> 'other'
        counts = get_base_counts("ATCGRYSWKM")
        assert counts['A'] == 1
        assert counts['T'] == 1
        assert counts['C'] == 1
        assert counts['G'] == 1
        assert counts['N'] == 0
        assert counts['other'] == 6
        assert counts['total'] == 10
    
    def test_empty_sequence(self):
        counts = get_base_counts("")
        assert counts['total'] == 0
        assert counts['A'] == 0
        assert counts['T'] == 0
        assert counts['C'] == 0
        assert counts['G'] == 0
        assert counts['N'] == 0
        assert counts['other'] == 0
    
    def test_whitespace_ignored(self):
        counts = get_base_counts("ATCG ATCG\nATCG")
        assert counts['A'] == 3
        assert counts['T'] == 3
        assert counts['C'] == 3
        assert counts['G'] == 3
        assert counts['total'] == 12


class TestGetSummaryStatistics:
    """Tests for get_summary_statistics function."""
    
    def test_comprehensive_stats(self):
        seq = "ATCGNNNATCG"
        stats = get_summary_statistics(seq)
        
        # Verify all keys are present
        assert 'length' in stats
        assert 'gc_content' in stats
        assert 'n_count' in stats
        assert 'at_gc_ratio' in stats
        assert 'base_counts' in stats
        
        # Verify values match individual function calls
        assert stats['length'] == length(seq)
        assert stats['gc_content'] == pytest.approx(gc_content(seq))
        assert stats['n_count'] == n_count(seq)
        assert stats['at_gc_ratio'] == pytest.approx(at_gc_ratio(seq))
        assert stats['base_counts'] == get_base_counts(seq)
    
    def test_balanced_sequence(self):
        seq = "ATCGATCG"
        stats = get_summary_statistics(seq)
        
        assert stats['length'] == 8
        assert stats['gc_content'] == pytest.approx(50.0)
        assert stats['n_count'] == 0
        assert stats['at_gc_ratio'] == pytest.approx(1.0)
        assert stats['base_counts']['total'] == 8
    
    def test_gc_rich_sequence(self):
        seq = "GGGGCCCC"
        stats = get_summary_statistics(seq)
        
        assert stats['length'] == 8
        assert stats['gc_content'] == pytest.approx(100.0)
        assert stats['at_gc_ratio'] == 0.0
    
    def test_at_rich_sequence(self):
        seq = "AAAAATTTTT"
        stats = get_summary_statistics(seq)
        
        assert stats['length'] == 10
        assert stats['gc_content'] == pytest.approx(0.0)
        assert stats['at_gc_ratio'] == float('inf')
    
    def test_rna_sequence(self):
        seq = "AUCGAUCG"
        stats = get_summary_statistics(seq)
        
        assert stats['length'] == 8
        assert stats['gc_content'] == pytest.approx(50.0)
        assert stats['base_counts']['T'] == 2  # U counted as T
