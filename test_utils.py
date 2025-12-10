#!/usr/bin/env python3
"""Pytest test suite for utils.py"""

import pytest
from pathlib import Path
from utils import (
    clean_sequence,
    safe_div,
    format_large_number,
    truncate_string,
    validate_file_path,
    calculate_percentage,
    is_valid_sequence,
    reverse_complement,
    get_file_size_mb,
    count_sequences_quick
)


class TestCleanSequence:
    """Tests for clean_sequence function."""
    
    def test_whitespace_and_newlines(self):
        assert clean_sequence("ATCG ATCG\nATCG") == "ATCGATCGATCG"
    
    def test_lowercase_to_uppercase(self):
        assert clean_sequence("  mklt vvga  ") == "MKLTVVGA"
    
    def test_spaces_between_chars(self):
        assert clean_sequence("A T C G") == "ATCG"
    
    def test_tabs_and_mixed_whitespace(self):
        assert clean_sequence("AT\tCG \n AT\tCG") == "ATCGATCG"
    
    def test_empty_string(self):
        assert clean_sequence("") == ""
    
    def test_only_whitespace(self):
        assert clean_sequence("   \n\t  ") == ""


class TestSafeDiv:
    """Tests for safe_div function."""
    
    def test_normal_division(self):
        assert safe_div(10, 2) == pytest.approx(5.0)
    
    def test_division_by_zero(self):
        assert safe_div(10, 0) == 0.0
    
    def test_zero_divided_by_nonzero(self):
        assert safe_div(0, 5) == 0.0
    
    def test_float_division(self):
        assert safe_div(7, 3) == pytest.approx(2.3333, rel=0.01)
    
    def test_both_zero(self):
        assert safe_div(0, 0) == 0.0


class TestFormatLargeNumber:
    """Tests for format_large_number function."""
    
    def test_integer_formatting(self):
        assert format_large_number(1234567) == "1,234,567"
    
    def test_float_with_default_decimals(self):
        assert format_large_number(1234.5678) == "1,234.57"
    
    def test_float_with_custom_decimals(self):
        assert format_large_number(999.999, decimals=0) == "1,000"
    
    def test_small_integer(self):
        assert format_large_number(42) == "42"
    
    def test_float_with_four_decimals(self):
        result = format_large_number(1234.56789, decimals=4)
        assert result == "1,234.5679"


class TestTruncateString:
    """Tests for truncate_string function."""
    
    def test_long_string_truncation(self):
        long_str = "This is a very long string that definitely needs to be truncated"
        result = truncate_string(long_str, 30)
        assert len(result) == 30
        assert result.endswith("...")
    
    def test_short_string_unchanged(self):
        short_str = "Short"
        result = truncate_string(short_str, 50)
        assert result == short_str
    
    def test_exact_length_unchanged(self):
        text = "X" * 50
        result = truncate_string(text, 50)
        assert result == text
    
    def test_custom_suffix(self):
        long_str = "A" * 100
        result = truncate_string(long_str, 20, suffix="[...]")
        assert len(result) == 20
        assert result.endswith("[...]")


class TestCalculatePercentage:
    """Tests for calculate_percentage function."""
    
    def test_simple_percentage(self):
        assert calculate_percentage(25, 100) == pytest.approx(25.0)
    
    def test_fractional_percentage(self):
        assert calculate_percentage(1, 3) == pytest.approx(33.33)
    
    def test_division_by_zero(self):
        assert calculate_percentage(10, 0) == 0.0
    
    def test_custom_decimals(self):
        result = calculate_percentage(2, 3, decimals=4)
        assert result == pytest.approx(66.6667)
    
    def test_zero_part(self):
        assert calculate_percentage(0, 100) == 0.0
    
    def test_100_percent(self):
        assert calculate_percentage(100, 100) == pytest.approx(100.0)


class TestIsValidSequence:
    """Tests for is_valid_sequence function."""
    
    def test_valid_dna_sequence(self):
        dna_alphabet = {'A', 'T', 'C', 'G', 'N'}
        assert is_valid_sequence("ATCGNATCG", dna_alphabet) is True
    
    def test_valid_case_insensitive(self):
        dna_alphabet = {'A', 'T', 'C', 'G', 'N'}
        assert is_valid_sequence("atcgnatcg", dna_alphabet) is True
    
    def test_invalid_character(self):
        dna_alphabet = {'A', 'T', 'C', 'G', 'N'}
        assert is_valid_sequence("ATCGX", dna_alphabet) is False
    
    def test_empty_sequence(self):
        dna_alphabet = {'A', 'T', 'C', 'G'}
        assert is_valid_sequence("", dna_alphabet) is False
    
    def test_whitespace_only(self):
        dna_alphabet = {'A', 'T', 'C', 'G'}
        assert is_valid_sequence("   \n\t  ", dna_alphabet) is False
    
    def test_protein_alphabet(self):
        protein_alphabet = set("ACDEFGHIKLMNPQRSTVWY")
        assert is_valid_sequence("MKLVVVGA", protein_alphabet) is True
        assert is_valid_sequence("MKLVVVGAX", protein_alphabet) is False


class TestReverseComplement:
    """Tests for reverse_complement function."""
    
    def test_simple_sequence(self):
        assert reverse_complement("ATCG") == "CGAT"
    
    def test_palindrome(self):
        assert reverse_complement("CGCG") == "CGCG"
    
    def test_all_a(self):
        assert reverse_complement("AAAAA") == "TTTTT"
    
    def test_complex_sequence(self):
        assert reverse_complement("GATTACA") == "TGTAATC"
    
    def test_rna_to_dna(self):
        # RNA input with U should give DNA output
        assert reverse_complement("AUCG") == "CGAT"
    
    def test_with_n(self):
        # N should map to N
        assert reverse_complement("ATCGN") == "NCGAT"
    
    def test_case_insensitive(self):
        result1 = reverse_complement("atcg")
        result2 = reverse_complement("ATCG")
        assert result1 == result2 == "CGAT"


class TestFileFunctions:
    """Tests for file-related functions."""
    
    def test_validate_file_path_existing(self):
        # test_mixed.fasta should exist in the project root
        assert validate_file_path("test_mixed.fasta") is True
    
    def test_validate_file_path_nonexistent(self):
        assert validate_file_path("nonexistent_file_xyz.fasta") is False
    
    def test_get_file_size_mb_existing(self):
        # test_mixed.fasta should have a size > 0
        size = get_file_size_mb("test_mixed.fasta")
        assert size > 0.0
    
    def test_get_file_size_mb_nonexistent(self):
        size = get_file_size_mb("nonexistent.txt")
        assert size == 0.0
    
    def test_count_sequences_quick_existing(self):
        # test_mixed.fasta has 12 sequences starting with '>'
        # (One sequence has a non-standard header detected by fallback logic)
        count = count_sequences_quick("test_mixed.fasta")
        assert count == 12
    
    def test_count_sequences_quick_nonexistent(self):
        count = count_sequences_quick("nonexistent.fasta")
        assert count == 0
    
    def test_count_sequences_quick_with_different_header_char(self):
        # With wrong header character, should find 0 sequences
        count = count_sequences_quick("test_mixed.fasta", header_char='@')
        assert count == 0


class TestFileFunctionsWithTmpPath:
    """Tests for file functions using temporary files."""
    
    def test_validate_file_path_with_tmp_file(self, tmp_path):
        # Create a temporary file
        temp_file = tmp_path / "test.txt"
        temp_file.write_text("content")
        
        assert validate_file_path(str(temp_file)) is True
    
    def test_get_file_size_with_known_content(self, tmp_path):
        # Create file with known content
        temp_file = tmp_path / "test.fasta"
        content = "A" * 1024  # 1KB
        temp_file.write_text(content)
        
        size_mb = get_file_size_mb(str(temp_file))
        assert size_mb == pytest.approx(0.001, rel=0.1)  # ~0.001 MB
    
    def test_count_sequences_quick_custom_fasta(self, tmp_path):
        # Create a simple FASTA file with 3 sequences
        temp_file = tmp_path / "test.fasta"
        content = ">seq1\nATCG\n>seq2\nGGCC\n>seq3\nAAAA\n"
        temp_file.write_text(content)
        
        count = count_sequences_quick(str(temp_file))
        assert count == 3
