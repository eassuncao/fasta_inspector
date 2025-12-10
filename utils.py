# AI Assistance: GitHub Copilot contributed utility function suggestions.

"""
Utility Functions Module

This module provides shared helper functions used across the bioinformatics project.
These utilities handle common tasks like sequence cleaning, safe mathematical operations,
and file handling.

Functions:
    clean_sequence(seq): Remove whitespace and convert to uppercase
    safe_div(numerator, denominator): Division with zero-denominator handling
    format_large_number(num): Format large numbers with comma separators
    truncate_string(s, max_len): Truncate strings with ellipsis
    validate_file_path(filepath): Check if a file exists and is readable
"""

import os
from typing import Union


def clean_sequence(seq: str) -> str:
    """
    Remove whitespace and newlines from a sequence and convert to uppercase.
    
    This is a common preprocessing step for biological sequences (DNA, RNA, protein)
    before analysis. Only alphabetic characters are retained.
    
    Args:
        seq (str): Biological sequence that may contain whitespace, newlines, or mixed case
        
    Returns:
        str: Cleaned sequence with only uppercase alphabetic characters
        
    Example:
        >>> clean_sequence('atcg ATCG\\nATCG')
        'ATCGATCGATCG'
        >>> clean_sequence('  M K L T  ')
        'MKLT'
    """
    # Filter out whitespace and convert to uppercase for consistent processing
    return ''.join(char.upper() for char in seq if char.isalpha())


def safe_div(numerator: Union[int, float], denominator: Union[int, float]) -> float:
    """
    Perform division with automatic handling of zero denominator.
    
    Returns 0.0 if the denominator is zero, avoiding ZeroDivisionError.
    Useful for calculating ratios and percentages where division by zero
    is a valid edge case that should return a sensible default.
    
    Args:
        numerator (int or float): The dividend
        denominator (int or float): The divisor
        
    Returns:
        float: Result of division, or 0.0 if denominator is zero
        
    Example:
        >>> safe_div(10, 2)
        5.0
        >>> safe_div(10, 0)
        0.0
        >>> safe_div(0, 0)
        0.0
    """
    if denominator == 0:
        return 0.0
    return numerator / denominator


def format_large_number(num: Union[int, float], decimals: int = 2) -> str:
    """
    Format large numbers with comma separators for readability.
    
    Integers are formatted with comma thousands separators. Floats are
    formatted with the specified number of decimal places.
    
    Args:
        num (int or float): Number to format
        decimals (int): Number of decimal places for floats (default: 2)
        
    Returns:
        str: Formatted number string with comma separators
        
    Example:
        >>> format_large_number(1234567)
        '1,234,567'
        >>> format_large_number(1234.5678, decimals=2)
        '1,234.57'
        >>> format_large_number(999.9, decimals=0)
        '1,000'
    """
    if isinstance(num, float):
        return f"{num:,.{decimals}f}"
    return f"{num:,}"


def truncate_string(s: str, max_len: int = 50, suffix: str = "...") -> str:
    """
    Truncate a string to a maximum length and add ellipsis.
    
    Useful for displaying long sequence headers or descriptions in reports.
    
    Args:
        s (str): String to truncate
        max_len (int): Maximum length including suffix (default: 50)
        suffix (str): Suffix to add when truncated (default: "...")
        
    Returns:
        str: Truncated string with suffix, or original if shorter than max_len
        
    Example:
        >>> truncate_string('This is a very long string that needs truncating', 20)
        'This is a very lo...'
        >>> truncate_string('Short', 20)
        'Short'
    """
    if len(s) <= max_len:
        return s
    return s[:max_len - len(suffix)] + suffix


def validate_file_path(filepath: str) -> bool:
    """
    Check if a file exists and is readable.
    
    Args:
        filepath (str): Path to the file to validate
        
    Returns:
        bool: True if file exists and is readable, False otherwise
        
    Example:
        >>> validate_file_path('existing_file.fasta')
        True
        >>> validate_file_path('nonexistent.txt')
        False
    """
    return os.path.isfile(filepath) and os.access(filepath, os.R_OK)


def calculate_percentage(part: Union[int, float], total: Union[int, float], 
                         decimals: int = 2) -> float:
    """
    Calculate percentage with safe division.
    
    Returns 0.0 if total is zero (same zero-division handling as safe_div).
    Useful for composition analysis and statistical calculations.
    
    Args:
        part (int or float): The part/subset value
        total (int or float): The total/whole value
        decimals (int): Number of decimal places to round to (default: 2)
        
    Returns:
        float: Percentage value (0.0 to 100.0), or 0.0 if total is zero
        
    Example:
        >>> calculate_percentage(25, 100)
        25.0
        >>> calculate_percentage(1, 3)
        33.33
        >>> calculate_percentage(10, 0)  # Safe division by zero
        0.0
        >>> calculate_percentage(2, 3, decimals=4)
        66.6667
    """
    if total == 0:
        return 0.0
    percentage = (part / total) * 100.0
    return round(percentage, decimals)


def is_valid_sequence(seq: str, alphabet: set) -> bool:
    """
    Check if a sequence contains only characters from a specified alphabet.
    
    Useful for validating sequences against a known set of valid characters
    (e.g., DNA alphabet, protein alphabet). The check is case-insensitive
    and ignores whitespace.
    
    Args:
        seq (str): Sequence to validate
        alphabet (set): Set of valid characters (case-insensitive)
        
    Returns:
        bool: True if all characters are in alphabet, False otherwise or if
              sequence is empty after cleaning
        
    Example:
        >>> is_valid_sequence('ATCG', {'A', 'T', 'C', 'G'})
        True
        >>> is_valid_sequence('ATCGX', {'A', 'T', 'C', 'G'})
        False
        >>> is_valid_sequence('atcg', {'A', 'T', 'C', 'G'})  # Case-insensitive
        True
    """
    cleaned = clean_sequence(seq)
    # Return False for empty sequences
    if not cleaned:
        return False
    return all(char in alphabet for char in cleaned)


def reverse_complement(seq: str) -> str:
    """
    Calculate the reverse complement of a DNA/RNA sequence.
    
    Note: U (uracil) in the input is treated as RNA, but the output reverse
    complement is always given in the DNA alphabet (A/T/C/G). Any U in the
    input is complemented to A (not U), resulting in DNA output.
    
    Args:
        seq (str): DNA or RNA sequence
        
    Returns:
        str: Reverse complement of the sequence (always DNA alphabet)
        
    Example:
        >>> reverse_complement('ATCG')
        'CGAT'
        >>> reverse_complement('AUCG')  # RNA input, DNA output
        'CGAT'
    """
    # Complement mapping for DNA and RNA
    complement = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'U': 'A',  # RNA
        'N': 'N',  # Ambiguous base
    }
    
    cleaned = clean_sequence(seq)
    
    # Build complement and reverse
    comp_seq = ''.join(complement.get(base, base) for base in cleaned)
    return comp_seq[::-1]


def get_file_size_mb(filepath: str) -> float:
    """
    Get file size in megabytes.
    
    Useful for checking file size before processing large files or
    displaying file information to users.
    
    Args:
        filepath (str): Path to the file
        
    Returns:
        float: File size in MB, or 0.0 if file doesn't exist or is not accessible
        
    Example:
        >>> get_file_size_mb('large_file.fasta')
        15.23
        >>> get_file_size_mb('nonexistent.txt')
        0.0
    """
    if not os.path.isfile(filepath):
        return 0.0
    try:
        size_bytes = os.path.getsize(filepath)
        return size_bytes / (1024 * 1024)
    except (OSError, IOError):
        # Return 0.0 if file size cannot be determined
        return 0.0


def count_sequences_quick(filepath: str, header_char: str = '>') -> int:
    """
    Quickly count sequences in a FASTA file by counting header lines.
    
    This is faster than parsing the entire file when you only need the count.
    Note: Only counts lines that start with header_char (default '>'). Will NOT
    detect sequences whose headers do not begin with this character.
    For files with non-standard headers, use fasta_reader.count_sequences() instead.
    
    Args:
        filepath (str): Path to the FASTA file
        header_char (str): Character that starts header lines (default: '>')
        
    Returns:
        int: Number of sequences (header lines) in the file, or 0 if file
             is not accessible or cannot be read
        
    Example:
        >>> count_sequences_quick('sequences.fasta')
        42
        >>> count_sequences_quick('empty.fasta')
        0
    """
    if not validate_file_path(filepath):
        return 0
    
    count = 0
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            for line in f:
                # Only count lines starting with header_char
                if line.strip().startswith(header_char):
                    count += 1
    except (IOError, UnicodeDecodeError):
        # Return 0 if file cannot be read
        return 0
    
    return count


if __name__ == "__main__":
    # Test/demonstration code
    print("Utility Functions Module - Test Suite\n")
    print("=" * 80)
    
    # Test clean_sequence
    print("\n1. clean_sequence()")
    print("-" * 80)
    test_seqs = [
        "ATCG ATCG\nATCG",
        "  mklt vvga  ",
        "A T C G",
    ]
    for seq in test_seqs:
        cleaned = clean_sequence(seq)
        print(f"Input: '{seq[:30]}...' → Output: '{cleaned}'")
    
    # Test safe_div
    print("\n2. safe_div()")
    print("-" * 80)
    divisions = [(10, 2), (10, 0), (0, 5), (7, 3)]
    for num, den in divisions:
        result = safe_div(num, den)
        print(f"{num} / {den} = {result:.2f}")
    
    # Test format_large_number
    print("\n3. format_large_number()")
    print("-" * 80)
    numbers = [1234567, 1234.5678, 999999.99]
    for num in numbers:
        formatted = format_large_number(num)
        print(f"{num} → {formatted}")
    
    # Test truncate_string
    print("\n4. truncate_string()")
    print("-" * 80)
    long_str = "This is a very long string that definitely needs to be truncated"
    print(f"Original ({len(long_str)} chars): {long_str}")
    print(f"Truncated (30): {truncate_string(long_str, 30)}")
    print(f"Truncated (50): {truncate_string(long_str, 50)}")
    
    # Test calculate_percentage
    print("\n5. calculate_percentage()")
    print("-" * 80)
    percentages = [(25, 100), (1, 3), (10, 0), (33, 99)]
    for part, total in percentages:
        pct = calculate_percentage(part, total)
        print(f"{part}/{total} = {pct}%")
    
    # Test is_valid_sequence
    print("\n6. is_valid_sequence()")
    print("-" * 80)
    dna_alphabet = {'A', 'T', 'C', 'G', 'N'}
    test_validation = [
        ("ATCGNATCG", dna_alphabet),
        ("ATCGX", dna_alphabet),
        ("ATCG", dna_alphabet),
    ]
    for seq, alphabet in test_validation:
        valid = is_valid_sequence(seq, alphabet)
        print(f"'{seq}' with DNA alphabet: {valid}")
    
    # Test reverse_complement
    print("\n7. reverse_complement()")
    print("-" * 80)
    rc_tests = ["ATCG", "AUCG", "GATTACA"]
    for seq in rc_tests:
        rc = reverse_complement(seq)
        print(f"{seq} → {rc}")
    
    print("\n" + "=" * 80)
    print("\nTo use this module, import the functions:")
    print("  from utils import clean_sequence, safe_div, truncate_string")
