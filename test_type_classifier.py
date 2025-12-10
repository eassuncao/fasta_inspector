#!/usr/bin/env python3
"""Pytest test suite for type_classifier.py"""

import pytest
from type_classifier import (
    SequenceClassificationError,
    detect_sequence_type,
    get_sequence_statistics,
    classify_fasta_sequences
)


class TestDetectSequenceType:
    """Tests for detect_sequence_type function."""
    
    def test_pure_dna(self):
        assert detect_sequence_type("ATCGATCGNNNN") == "DNA"
    
    def test_pure_protein(self):
        assert detect_sequence_type("MKLVVVGAGGKS") == "PROTEIN"
    
    def test_dna_with_high_percentage(self):
        # ATCGATCGATCG - all DNA core alphabet
        assert detect_sequence_type("ATCGATCGATCG") == "DNA"
    
    def test_protein_with_overlap_chars(self):
        # MKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEY
        # Should be classified as PROTEIN despite some overlap
        seq = "MKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEY"
        assert detect_sequence_type(seq) == "PROTEIN"
    
    def test_mixed_below_threshold(self):
        # ATCGMKLVVV: 4 DNA (40%), 6 protein (60%) -> PROTEIN
        assert detect_sequence_type("ATCGMKLVVV") == "PROTEIN"
    
    def test_rna_with_uracil(self):
        # AUCGAUCGAUC - RNA should be classified as DNA
        assert detect_sequence_type("AUCGAUCGAUC") == "DNA"
    
    def test_case_insensitive(self):
        assert detect_sequence_type("atcgatcg") == "DNA"
        assert detect_sequence_type("mklvvvga") == "PROTEIN"
    
    def test_with_whitespace(self):
        assert detect_sequence_type("ATCG ATCG\nATCG") == "DNA"
    
    def test_empty_sequence_raises(self):
        with pytest.raises(SequenceClassificationError, match="Empty or invalid"):
            detect_sequence_type("")
    
    def test_none_input_raises(self):
        with pytest.raises(SequenceClassificationError, match="Empty or invalid"):
            detect_sequence_type(None)
    
    def test_non_string_input_raises(self):
        with pytest.raises(SequenceClassificationError, match="Empty or invalid"):
            detect_sequence_type(123)
    
    def test_only_whitespace_raises(self):
        with pytest.raises(SequenceClassificationError, match="no valid alphabetic"):
            detect_sequence_type("   \n\t  ")


class TestGetSequenceStatistics:
    """Tests for get_sequence_statistics function."""
    
    def test_dna_statistics(self):
        seq = "ATCGATCG"
        stats = get_sequence_statistics(seq)
        
        assert 'classification' in stats
        assert stats['classification'] == "DNA"
        assert 'total_chars' in stats
        assert 'valid_chars' in stats
        assert 'dna_chars' in stats
        assert 'dna_proportion' in stats
        
        # Check that dna_proportion is in valid range
        assert 0.0 <= stats['dna_proportion'] <= 1.0
    
    def test_protein_statistics(self):
        seq = "MKLVVVGAGGKS"
        stats = get_sequence_statistics(seq)
        
        assert stats['classification'] == "PROTEIN"
        assert 0.0 <= stats['dna_proportion'] <= 1.0
    
    def test_statistics_length(self):
        seq = "ATCGATCG"
        stats = get_sequence_statistics(seq)
        
        assert stats['valid_chars'] == 8
        assert stats['total_chars'] == 8
    
    def test_statistics_with_whitespace(self):
        seq = "ATCG ATCG"
        stats = get_sequence_statistics(seq)
        
        assert stats['total_chars'] == 9  # includes space
        assert stats['valid_chars'] == 8  # alphabetic only
    
    def test_dna_proportion_high(self):
        seq = "ATCGATCGNNNN"
        stats = get_sequence_statistics(seq)
        
        # All characters are DNA core alphabet
        assert stats['dna_proportion'] == pytest.approx(1.0)
    
    def test_dna_proportion_low(self):
        seq = "MKLVVVGAGGKS"
        stats = get_sequence_statistics(seq)
        
        # Protein sequence should have low DNA proportion
        assert stats['dna_proportion'] < 0.9
    
    def test_empty_sequence_raises(self):
        with pytest.raises(SequenceClassificationError):
            get_sequence_statistics("")
    
    def test_none_input_raises(self):
        with pytest.raises(SequenceClassificationError):
            get_sequence_statistics(None)


class TestClassifyFastaSequences:
    """Tests for classify_fasta_sequences function."""
    
    def test_classify_mixed_sequences(self):
        # Prepare test data
        fasta_data = [
            ("seq1", "ATCGATCGATCG"),
            ("seq2", "MKLVVVGAGGKS"),
            ("seq3", "GGCCGGCCGGCC"),
            ("seq4", "ACDEFGHIKLMN")
        ]
        
        results = classify_fasta_sequences(fasta_data)
        
        assert len(results) == 4
        
        # Check structure of results
        for header, sequence, classification in results:
            assert isinstance(header, str)
            assert isinstance(sequence, str)
            assert classification in ["DNA", "PROTEIN"] or classification.startswith("ERROR:")
    
    def test_classify_all_dna(self):
        fasta_data = [
            ("dna1", "ATCGATCG"),
            ("dna2", "GGCCGGCC"),
            ("dna3", "AAAATTTT")
        ]
        
        results = classify_fasta_sequences(fasta_data)
        
        for _, _, classification in results:
            assert classification == "DNA"
    
    def test_classify_all_protein(self):
        fasta_data = [
            ("prot1", "MKLVVVGA"),
            ("prot2", "ACDEFGHI"),
            ("prot3", "PQRSTVWY")
        ]
        
        results = classify_fasta_sequences(fasta_data)
        
        for _, _, classification in results:
            assert classification == "PROTEIN"
    
    def test_classify_with_error(self):
        # Include an empty sequence that will cause an error
        fasta_data = [
            ("seq1", "ATCGATCG"),
            ("empty", ""),
            ("seq3", "MKLVVVGA")
        ]
        
        results = classify_fasta_sequences(fasta_data)
        
        assert len(results) == 3
        # Second result should be an error
        assert results[1][2].startswith("ERROR:")
    
    def test_empty_list(self):
        results = classify_fasta_sequences([])
        assert results == []


class TestIntegrationWithRealFile:
    """Integration tests using test_mixed.fasta"""
    
    def test_classify_test_mixed_fasta(self):
        from fasta_reader import read_fasta
        
        # Read sequences from test file
        sequences = list(read_fasta("test_mixed.fasta"))
        
        # Classify all sequences
        results = classify_fasta_sequences(sequences)
        
        # Should have 13 sequences
        assert len(results) == 13
        
        # Count DNA vs PROTEIN
        dna_count = sum(1 for _, _, t in results if t == "DNA")
        protein_count = sum(1 for _, _, t in results if t == "PROTEIN")
        
        # Based on current behavior: 4 DNA, 9 protein
        assert dna_count == 4
        assert protein_count == 9
    
    def test_individual_classification(self):
        from fasta_reader import read_fasta
        
        sequences = list(read_fasta("test_mixed.fasta"))
        
        # Test that each sequence gets classified
        for header, sequence in sequences:
            seq_type = detect_sequence_type(sequence)
            assert seq_type in ["DNA", "PROTEIN"]
            
            # Also test statistics
            stats = get_sequence_statistics(sequence)
            assert stats['classification'] in ["DNA", "PROTEIN"]
