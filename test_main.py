#!/usr/bin/env python3
"""Pytest test suite for main.py"""

import pytest
from main import (
    analyze_fasta_file,
    print_summary
)


class TestAnalyzeFastaFile:
    """Tests for analyze_fasta_file function."""
    
    def test_analyze_test_mixed_fasta(self):
        total, dna, protein, total_bases, total_residues = analyze_fasta_file("test_mixed.fasta")
        
        # Based on current behavior
        assert total == 13
        assert dna == 4
        assert protein == 9
        assert total_bases > 0
        assert total_residues > 0
    
    def test_analyze_returns_tuple(self):
        result = analyze_fasta_file("test_mixed.fasta")
        
        assert isinstance(result, tuple)
        assert len(result) == 5
        
        total, dna, protein, total_bases, total_residues = result
        assert isinstance(total, int)
        assert isinstance(dna, int)
        assert isinstance(protein, int)
        assert isinstance(total_bases, int)
        assert isinstance(total_residues, int)
    
    def test_analyze_nonexistent_file(self):
        # Should handle error gracefully and return zeros
        total, dna, protein, total_bases, total_residues = analyze_fasta_file("nonexistent_xyz.fasta")
        
        assert total == 0
        assert dna == 0
        assert protein == 0
        assert total_bases == 0
        assert total_residues == 0
    
    def test_analyze_counts_add_up(self):
        total, dna, protein, total_bases, total_residues = analyze_fasta_file("test_mixed.fasta")
        
        # DNA + protein should equal or be less than total
        # (could be less if there are unclassified sequences)
        assert dna + protein <= total


class TestPrintSummary:
    """Tests for print_summary function using capsys."""
    
    def test_print_summary_output(self, capsys):
        print_summary(13, 4, 9)
        
        captured = capsys.readouterr()
        output = captured.out
        
        # Check for key elements in output
        assert "Total Sequences:" in output
        assert "13" in output
        assert "DNA/RNA Sequences:" in output
        assert "4" in output
        assert "Protein Sequences:" in output
        assert "9" in output
    
    def test_print_summary_percentages(self, capsys):
        print_summary(13, 4, 9)
        
        captured = capsys.readouterr()
        output = captured.out
        
        # Check for percentages
        assert "30.8%" in output  # 4/13 * 100
        assert "69.2%" in output  # 9/13 * 100
    
    def test_print_summary_with_zeros(self, capsys):
        print_summary(0, 0, 0)
        
        captured = capsys.readouterr()
        output = captured.out
        
        # Should handle zero total gracefully
        assert "Total Sequences:" in output
        assert "0" in output
    
    def test_print_summary_with_unclassified(self, capsys):
        # 10 total, 3 DNA, 5 protein = 2 unclassified
        print_summary(10, 3, 5)
        
        captured = capsys.readouterr()
        output = captured.out
        
        # Should show unclassified count
        assert "Unclassified:" in output
        assert "2" in output


class TestIntegrationWithRealFile:
    """Integration tests with actual FASTA file."""
    
    def test_full_analysis_workflow(self):
        # Run the full analysis
        total, dna, protein, total_bases, total_residues = analyze_fasta_file("test_mixed.fasta")
        
        # Verify counts
        assert total > 0
        assert dna > 0
        assert protein > 0
        assert dna + protein == total  # No unclassified in test_mixed.fasta
        assert total_bases > 0
        assert total_residues > 0
    
    def test_analysis_consistent_with_file(self):
        from fasta_reader import count_sequences
        
        # Count using fasta_reader
        file_count = count_sequences("test_mixed.fasta")
        
        # Count using analyze_fasta_file
        total, _, _, _, _ = analyze_fasta_file("test_mixed.fasta")
        
        # Should match
        assert total == file_count
    
    def test_analysis_matches_classification(self):
        from fasta_reader import read_fasta
        from type_classifier import detect_sequence_type
        
        # Manually classify sequences
        sequences = list(read_fasta("test_mixed.fasta"))
        manual_dna = 0
        manual_protein = 0
        
        for header, sequence in sequences:
            seq_type = detect_sequence_type(sequence)
            if seq_type == "DNA":
                manual_dna += 1
            elif seq_type == "PROTEIN":
                manual_protein += 1
        
        # Compare with analyze_fasta_file results
        total, dna, protein, total_bases, total_residues = analyze_fasta_file("test_mixed.fasta")
        
        assert dna == manual_dna
        assert protein == manual_protein
