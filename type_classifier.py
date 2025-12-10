# AI Assistance: GitHub Copilot helped generate this classification module.

"""
Sequence Type Classification Module

This module provides functionality for classifying biological sequences as
either DNA or protein based on their nucleotide composition.

Functions:
    detect_sequence_type(seq): Classify a sequence as 'DNA' or 'PROTEIN'
"""

from typing import Set, List, Tuple


# Core DNA alphabet for classification
# Uses only standard bases plus N (any base) and U (RNA uracil)
# This avoids false positives from amino acids that overlap with IUPAC ambiguity codes
DNA_CORE_ALPHABET: Set[str] = set('ACGTUN')

# Full IUPAC nucleotide alphabet (for reference/validation if needed)
# Ambiguity codes: R (purine), Y (pyrimidine), S (strong), W (weak),
#                  K (keto), M (amino), B (not A), D (not C), H (not G), V (not T)
DNA_IUPAC_ALPHABET: Set[str] = set('ACGTUNRYSWKMBDHV')

# Threshold for DNA classification (90% of valid characters must be core DNA)
DNA_THRESHOLD: float = 0.90


class SequenceClassificationError(Exception):
    """Exception raised when a sequence cannot be classified."""
    pass


def detect_sequence_type(seq: str) -> str:
    """
    Classify a biological sequence as DNA or protein.
    
    This function analyzes the composition of a sequence and determines whether
    it is DNA or protein based on the proportion of DNA-specific characters.
    
    Classification logic:
    - If ≥90% of valid characters belong to the core DNA alphabet (ACGTUN),
      the sequence is classified as DNA.
    - Otherwise, it is classified as protein.
    - Uses core DNA alphabet to avoid false positives from amino acids that
      overlap with IUPAC ambiguity codes (e.g., M, K, V, etc.)
    - Whitespace and newline characters are automatically ignored.
    - The comparison is case-insensitive.
    
    Args:
        seq (str): The biological sequence to classify
        
    Returns:
        str: Either 'DNA' or 'PROTEIN'
        
    Raises:
        SequenceClassificationError: If the sequence is empty, contains only
                                    whitespace, or has no valid alphabetic
                                    characters
        
    Example:
        >>> detect_sequence_type('ATCGATCG')
        'DNA'
        >>> detect_sequence_type('MKLVVVGAGGVGKSALT')
        'PROTEIN'
        >>> detect_sequence_type('atcg')  # Case-insensitive
        'DNA'
    """
    # Validate input
    if not seq or not isinstance(seq, str):
        raise SequenceClassificationError("Empty or invalid sequence provided")
    
    # Convert to uppercase for case-insensitive comparison
    seq_upper = seq.upper()
    
    # Filter out whitespace and newline characters
    # Keep only alphabetic characters for analysis
    valid_chars = [char for char in seq_upper if char.isalpha()]
    
    # Check if we have any valid characters after filtering
    if not valid_chars:
        raise SequenceClassificationError(
            "Sequence contains no valid alphabetic characters"
        )
    
    total_valid = len(valid_chars)
    
    # Count how many characters belong to the core DNA alphabet
    # Using core alphabet (ACGTUN) avoids false positives from overlapping amino acids
    dna_char_count = sum(1 for char in valid_chars if char in DNA_CORE_ALPHABET)
    
    # Calculate the proportion of DNA characters
    dna_proportion = dna_char_count / total_valid
    
    # Apply the 90% threshold rule
    if dna_proportion >= DNA_THRESHOLD:
        return 'DNA'
    else:
        return 'PROTEIN'


def get_sequence_statistics(seq: str) -> dict:
    """
    Get detailed statistics about a sequence's composition.
    
    This helper function provides insights into what characters are present
    in the sequence and their proportions, useful for understanding why
    a particular classification was made.
    
    Args:
        seq (str): The biological sequence to analyze
        
    Returns:
        dict: A dictionary containing:
            - 'total_chars': Total number of characters (including whitespace)
            - 'valid_chars': Number of alphabetic characters
            - 'dna_chars': Number of DNA alphabet characters
            - 'dna_proportion': Proportion of DNA characters (0.0 to 1.0)
            - 'classification': The detected sequence type ('DNA' or 'PROTEIN')
            
    Raises:
        SequenceClassificationError: If the sequence has no valid characters
        
    Example:
        >>> stats = get_sequence_statistics('ATCGATCG')
        >>> stats['classification']
        'DNA'
        >>> stats['dna_proportion']
        1.0
    """
    # Validate input
    if not seq or not isinstance(seq, str):
        raise SequenceClassificationError("Empty or invalid sequence provided")
    
    seq_upper = seq.upper()
    valid_chars = [char for char in seq_upper if char.isalpha()]
    
    if not valid_chars:
        raise SequenceClassificationError(
            "Sequence contains no valid alphabetic characters"
        )
    
    total_valid = len(valid_chars)
    dna_char_count = sum(1 for char in valid_chars if char in DNA_CORE_ALPHABET)
    dna_proportion = dna_char_count / total_valid
    
    return {
        'total_chars': len(seq),
        'valid_chars': total_valid,
        'dna_chars': dna_char_count,
        'dna_proportion': dna_proportion,
        'classification': 'DNA' if dna_proportion >= DNA_THRESHOLD else 'PROTEIN'
    }


def classify_fasta_sequences(
    fasta_data: List[Tuple[str, str]]
) -> List[Tuple[str, str, str]]:
    """
    Classify multiple FASTA sequences in batch.
    
    This function processes a list of (header, sequence) tuples and returns
    the classification for each sequence. If classification fails for any
    sequence, the error message is included in the result.
    
    Args:
        fasta_data (List[Tuple[str, str]]): List of (header, sequence) tuples
        
    Returns:
        List[Tuple[str, str, str]]: List of (header, sequence, type) tuples
                                     where type is 'DNA', 'PROTEIN', or
                                     'ERROR: <message>' for failed classifications
        
    Example:
        >>> data = [('seq1', 'ATCG'), ('seq2', 'MKLV')]
        >>> results = classify_fasta_sequences(data)
        >>> len(results)
        2
    """
    results = []
    for header, sequence in fasta_data:
        try:
            seq_type = detect_sequence_type(sequence)
            results.append((header, sequence, seq_type))
        except SequenceClassificationError as e:
            # Include error message in result for debugging
            results.append((header, sequence, f"ERROR: {str(e)}"))
    return results


if __name__ == "__main__":
    # Test/demonstration code
    import sys
    
    # Try to import fasta_reader if available
    try:
        from fasta_reader import read_fasta
        HAS_FASTA_READER = True
    except ImportError:
        HAS_FASTA_READER = False
        print("Warning: fasta_reader module not found. Limited functionality.\n")
    
    # If a FASTA file is provided, process it
    if len(sys.argv) > 1 and HAS_FASTA_READER:
        fasta_file = sys.argv[1]
        try:
            print(f"Classifying sequences from: {fasta_file}\n")
            print("=" * 80)
            
            dna_count = 0
            protein_count = 0
            
            for i, (header, sequence) in enumerate(read_fasta(fasta_file), 1):
                seq_type = detect_sequence_type(sequence)
                stats = get_sequence_statistics(sequence)
                
                print(f"\nSequence {i}: {seq_type}")
                print(f"Header: {header[:60]}{'...' if len(header) > 60 else ''}")
                print(f"Length: {len(sequence)} characters")
                print(f"DNA proportion: {stats['dna_proportion']:.2%}")
                print(f"Preview: {sequence[:50]}{'...' if len(sequence) > 50 else ''}")
                
                if seq_type == 'DNA':
                    dna_count += 1
                else:
                    protein_count += 1
            
            print("\n" + "=" * 80)
            print(f"\nSummary:")
            print(f"  DNA sequences: {dna_count}")
            print(f"  Protein sequences: {protein_count}")
            print(f"  Total: {dna_count + protein_count}")
            
        except Exception as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)
    
    # Otherwise, run basic tests
    else:
        print("Type Classifier Module - Basic Tests\n")
        print("=" * 80)
        
        test_cases = [
            ("Pure DNA", "ATCGATCGATCG"),
            ("DNA with ambiguity (N)", "ATCGNNNATCG"),
            ("DNA lowercase", "atcgatcg"),
            ("RNA with uracil", "AUCGAUCGAUC"),
            ("Pure protein", "MKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEY"),
            ("Protein with common amino acids", "MGSSHHHHHHSSGLVPRGSHMLE"),
            ("Edge case - DNA with one protein AA", "ATCGATCGATCGM"),  # 92.3% DNA (12/13)
            ("Edge case - mixed below threshold", "ATCGMKLVVV"),  # 40% DNA (4/10)
            ("Edge case - exactly at threshold", "ATCGATCGATL"),  # 90.9% DNA (10/11)
        ]
        
        for name, sequence in test_cases:
            try:
                result = detect_sequence_type(sequence)
                stats = get_sequence_statistics(sequence)
                print(f"\nTest: {name}")
                print(f"Sequence: {sequence[:50]}{'...' if len(sequence) > 50 else ''}")
                print(f"Classification: {result}")
                print(f"DNA proportion: {stats['dna_proportion']:.2%}")
            except SequenceClassificationError as e:
                print(f"\nTest: {name}")
                print(f"Error: {e}")
        
        print("\n" + "=" * 80)
        
        if HAS_FASTA_READER:
            print("\nTo test with a FASTA file, run:")
            print("  python type_classifier.py <fasta_file>")
            print("\nExample:")
            print("  python type_classifier.py test_mixed.fasta")
        else:
            print("\nNote: Install fasta_reader.py in the same directory to test with FASTA files.")
