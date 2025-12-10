# AI Assistance: GitHub Copilot provided partial code suggestions for DNA metric functions.

"""
Nucleic Acid Sequence Metrics Module

This module provides functions for calculating various metrics on nucleic acid sequences
(DNA or RNA). All functions assume the input sequence has already been validated.

Functions:
    gc_content(seq): Calculate GC content percentage
    n_count(seq): Count ambiguous 'N' bases
    at_gc_ratio(seq): Calculate AT/GC ratio
    length(seq): Get sequence length (excluding whitespace)
"""


def _clean_sequence(seq: str) -> str:
    """
    Remove whitespace and newlines from a sequence.
    
    Helper function used by all metric functions to normalize input.
    
    Args:
        seq (str): DNA/RNA sequence that may contain whitespace or newlines
        
    Returns:
        str: Cleaned sequence with only alphabetic characters (uppercase)
    """
    # Filter out whitespace and convert to uppercase for consistent processing
    return ''.join(char.upper() for char in seq if char.isalpha())


def length(seq: str) -> int:
    """
    Calculate the length of a DNA/RNA sequence.
    
    Whitespace and newlines are automatically removed before counting.
    
    Args:
        seq (str): DNA/RNA sequence
        
    Returns:
        int: Number of bases in the sequence (excluding whitespace)
        
    Example:
        >>> length('ATCG ATCG')
        8
        >>> length('ATCG\\nATCG')
        8
    """
    cleaned = _clean_sequence(seq)
    return len(cleaned)


def gc_content(seq: str) -> float:
    """
    Calculate the GC content of a DNA/RNA sequence.
    
    GC content is the percentage of bases that are either Guanine (G) or 
    Cytosine (C). This is an important metric as GC-rich regions have 
    different properties than AT-rich regions.
    Note: U (uracil) is included in the total base count but not in GC count.
    
    Args:
        seq (str): DNA/RNA sequence (whitespace is ignored)
        
    Returns:
        float: GC content as a percentage (0.0 to 100.0)
               Returns 0.0 if sequence is empty or contains no valid bases
        
    Example:
        >>> gc_content('ATCG')
        50.0
        >>> gc_content('AAAA')
        0.0
        >>> gc_content('GGCC')
        100.0
    """
    cleaned = _clean_sequence(seq)
    
    # Handle empty sequence
    if not cleaned:
        return 0.0
    
    # Count G and C bases
    gc_count = sum(1 for base in cleaned if base in 'GC')
    
    # Calculate percentage
    gc_percentage = (gc_count / len(cleaned)) * 100.0
    
    return gc_percentage


def n_count(seq: str) -> int:
    """
    Count the number of ambiguous 'N' bases in a DNA/RNA sequence.
    
    'N' represents any nucleotide (A, C, G, T, or U) and is used when the 
    base at a position is unknown or ambiguous. A high N count may indicate 
    low sequencing quality.
    
    Args:
        seq (str): DNA/RNA sequence (whitespace is ignored)
        
    Returns:
        int: Number of 'N' bases in the sequence
        
    Example:
        >>> n_count('ATCGNNNATCG')
        3
        >>> n_count('ATCGATCG')
        0
    """
    cleaned = _clean_sequence(seq)
    
    # Count occurrences of 'N'
    return cleaned.count('N')


def at_gc_ratio(seq: str) -> float:
    """
    Calculate the AT/GC ratio of a DNA/RNA sequence.
    
    This ratio compares the proportion of AT (or AU for RNA) bases to GC bases. 
    A ratio > 1 indicates more AT bases, while < 1 indicates more GC bases.
    Note: U (uracil) is treated as equivalent to T (thymine) for RNA sequences.
    
    Args:
        seq (str): DNA/RNA sequence (whitespace is ignored)
        
    Returns:
        float: AT/GC ratio. Returns:
               - 0.0 if sequence is empty
               - float('inf') if there are AT bases but no GC bases
               - 0.0 if there are GC bases but no AT bases
               - The calculated ratio otherwise
        
    Example:
        >>> at_gc_ratio('AAAATTTT')
        inf
        >>> at_gc_ratio('GGGGCCCC')
        0.0
        >>> at_gc_ratio('ATCG')
        1.0
        >>> at_gc_ratio('AAATCG')
        1.5
        >>> at_gc_ratio('AUCG')  # RNA with uracil
        1.0
    """
    cleaned = _clean_sequence(seq)
    
    # Handle empty sequence
    if not cleaned:
        return 0.0
    
    # Count AT and GC bases (U is treated as T for RNA sequences)
    at_count = sum(1 for base in cleaned if base in 'ATU')
    gc_count = sum(1 for base in cleaned if base in 'GC')
    
    # Handle edge case: no GC bases (division by zero)
    if gc_count == 0:
        # If there are AT bases but no GC, return infinity
        if at_count > 0:
            return float('inf')
        # If there are neither AT nor GC (only N's or other ambiguity codes)
        else:
            return 0.0
    
    # Calculate ratio
    ratio = at_count / gc_count
    
    return ratio


def get_base_counts(seq: str) -> dict:
    """
    Get detailed counts of all bases in a DNA/RNA sequence.
    
    Helper function that provides a comprehensive breakdown of base composition.
    Note: U (uracil) is counted together with T (thymine) for consistency.
    
    Args:
        seq (str): DNA/RNA sequence (whitespace is ignored)
        
    Returns:
        dict: Dictionary with counts for each base type:
            - 'A': Adenine count
            - 'T': Thymine count (includes U/uracil for RNA)
            - 'C': Cytosine count
            - 'G': Guanine count
            - 'N': Ambiguous base count
            - 'other': Count of other IUPAC ambiguity codes
            - 'total': Total sequence length
            
    Example:
        >>> get_base_counts('ATCGNNNATCG')
        {'A': 2, 'T': 2, 'C': 2, 'G': 2, 'N': 3, 'other': 0, 'total': 11}
        >>> get_base_counts('AUCG')  # RNA with uracil
        {'A': 1, 'T': 1, 'C': 1, 'G': 1, 'N': 0, 'other': 0, 'total': 4}
    """
    cleaned = _clean_sequence(seq)
    
    # Initialize counts
    counts = {
        'A': 0,
        'T': 0,  # T includes U (uracil) for RNA sequences
        'C': 0,
        'G': 0,
        'N': 0,
        'other': 0,  # Other IUPAC ambiguity codes
        'total': len(cleaned)
    }
    
    # Count each base
    for base in cleaned:
        if base == 'U':
            # Treat U (uracil) as T (thymine) for RNA sequences
            counts['T'] += 1
        elif base in counts:
            counts[base] += 1
        else:
            # Other ambiguity codes (R, Y, S, W, K, M, B, D, H, V)
            counts['other'] += 1
    
    return counts


def get_summary_statistics(seq: str) -> dict:
    """
    Get comprehensive summary statistics for a DNA/RNA sequence.
    
    Convenience function that calculates all metrics at once.
    
    Args:
        seq (str): DNA/RNA sequence (whitespace is ignored)
        
    Returns:
        dict: Dictionary containing:
            - 'length': Sequence length
            - 'gc_content': GC content percentage
            - 'n_count': Number of N bases
            - 'at_gc_ratio': AT/GC ratio
            - 'base_counts': Detailed base composition
            
    Example:
        >>> stats = get_summary_statistics('ATCGNNNATCG')
        >>> stats['length']
        11
        >>> stats['gc_content']
        36.36...
    """
    return {
        'length': length(seq),
        'gc_content': gc_content(seq),
        'n_count': n_count(seq),
        'at_gc_ratio': at_gc_ratio(seq),
        'base_counts': get_base_counts(seq)
    }


if __name__ == "__main__":
    # Test/demonstration code
    print("DNA Metrics Module - Test Suite\n")
    print("=" * 80)
    
    test_sequences = [
        ("Pure DNA - balanced", "ATCGATCG"),
        ("GC-rich", "GGGGCCCCGGCC"),
        ("AT-rich", "AAAAATTTTT"),
        ("With ambiguity", "ATCGNNNATCG"),
        ("With whitespace", "ATCG ATCG\nATCG"),
        ("Only N's", "NNNNNN"),
        ("Empty", ""),
        ("Mixed ambiguity", "ATCGRYSWKM"),
    ]
    
    for name, seq in test_sequences:
        print(f"\nTest: {name}")
        print(f"Sequence: '{seq}'")
        
        stats = get_summary_statistics(seq)
        
        print(f"Length: {stats['length']}")
        print(f"GC Content: {stats['gc_content']:.2f}%")
        print(f"N Count: {stats['n_count']}")
        
        ratio = stats['at_gc_ratio']
        if ratio == float('inf'):
            print(f"AT/GC Ratio: inf (no GC bases)")
        else:
            print(f"AT/GC Ratio: {ratio:.2f}")
        
        print(f"Base Counts: A={stats['base_counts']['A']}, "
              f"T={stats['base_counts']['T']}, "
              f"C={stats['base_counts']['C']}, "
              f"G={stats['base_counts']['G']}, "
              f"N={stats['base_counts']['N']}, "
              f"Other={stats['base_counts']['other']}")
    
    print("\n" + "=" * 80)
    print("\nTo use this module, import the functions:")
    print("  from dna_metrics import gc_content, n_count, at_gc_ratio, length")
