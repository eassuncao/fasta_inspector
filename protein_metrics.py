# AI Assistance: GitHub Copilot assisted in generating protein metric functions.

"""
Protein Sequence Metrics Module

This module provides functions for calculating various metrics on protein sequences.
All functions assume the input sequence has already been validated as protein.

Functions:
    amino_acid_composition(seq): Count each amino acid in the sequence
    molecular_weight(seq): Calculate approximate molecular weight in Daltons
    length(seq): Get sequence length (excluding whitespace)
"""

from typing import Dict, Any


# Molecular weights of amino acids (average isotopic composition, in Daltons)
# Based on the monoisotopic mass of each residue in a polypeptide chain
# (i.e., after peptide bond formation, minus H2O)
AMINO_ACID_WEIGHTS: Dict[str, float] = {
    'A': 89.09,   # Alanine
    'R': 174.20,  # Arginine
    'N': 132.12,  # Asparagine
    'D': 133.10,  # Aspartic acid
    'C': 121.15,  # Cysteine
    'E': 147.13,  # Glutamic acid
    'Q': 146.15,  # Glutamine
    'G': 75.07,   # Glycine
    'H': 155.16,  # Histidine
    'I': 131.17,  # Isoleucine
    'L': 131.17,  # Leucine
    'K': 146.19,  # Lysine
    'M': 149.21,  # Methionine
    'F': 165.19,  # Phenylalanine
    'P': 115.13,  # Proline
    'S': 105.09,  # Serine
    'T': 119.12,  # Threonine
    'W': 204.23,  # Tryptophan
    'Y': 181.19,  # Tyrosine
    'V': 117.15,  # Valine
}

# Standard 20 amino acids (for validation)
STANDARD_AMINO_ACIDS = set(AMINO_ACID_WEIGHTS.keys())


class ProteinMetricsError(Exception):
    """Exception raised for errors in protein sequence metrics calculation."""
    pass


def _clean_sequence(seq: str) -> str:
    """
    Remove whitespace and newlines from a sequence.
    
    Helper function used by all metric functions to normalize input.
    Ensures consistent processing across all protein metrics calculations.
    
    Args:
        seq (str): Protein sequence that may contain whitespace or newlines
        
    Returns:
        str: Cleaned sequence with only alphabetic characters (uppercase)
        
    Example:
        >>> _clean_sequence('MKLT \\n VVGA')
        'MKLTVVGA'
    """
    # Filter out whitespace and convert to uppercase for consistent processing
    return ''.join(char.upper() for char in seq if char.isalpha())


def length(seq: str) -> int:
    """
    Calculate the length of a protein sequence.
    
    Whitespace and newlines are automatically removed before counting.
    
    Args:
        seq (str): Protein sequence
        
    Returns:
        int: Number of amino acid residues in the sequence (excluding whitespace)
        
    Raises:
        ProteinMetricsError: If sequence is empty or contains no valid characters
        
    Example:
        >>> length('MKLT VVGA')
        8
        >>> length('MKLT\\nVVGA')
        8
    """
    cleaned = _clean_sequence(seq)
    
    if not cleaned:
        raise ProteinMetricsError("Sequence is empty or contains no valid characters")
    
    return len(cleaned)


def amino_acid_composition(seq: str) -> Dict[str, int]:
    """
    Calculate the amino acid composition of a protein sequence.
    
    Returns a dictionary with counts for each of the 20 standard amino acids
    present in the sequence. Only uppercase amino acids are counted after
    cleaning the sequence. Non-standard amino acids are ignored.
    
    Args:
        seq (str): Protein sequence (whitespace is ignored)
        
    Returns:
        dict: Dictionary mapping amino acid one-letter codes to their counts.
              Only amino acids present in the sequence are included.
              Only the 20 standard amino acids (A,R,N,D,C,E,Q,G,H,I,L,K,M,F,P,S,T,W,Y,V) are counted.
              
    Raises:
        ProteinMetricsError: If sequence is empty or contains no valid amino acids
        
    Example:
        >>> amino_acid_composition('MKLVVVGA')
        {'M': 1, 'K': 1, 'L': 1, 'V': 3, 'G': 1, 'A': 1}
        >>> amino_acid_composition('AAA')
        {'A': 3}
    """
    cleaned = _clean_sequence(seq)
    
    if not cleaned:
        raise ProteinMetricsError("Sequence is empty or contains no valid characters")
    
    # Initialize composition dictionary (only for amino acids that appear)
    composition: Dict[str, int] = {}
    
    # Count each amino acid (only standard 20)
    for aa in cleaned:
        if aa in STANDARD_AMINO_ACIDS:
            composition[aa] = composition.get(aa, 0) + 1
    
    # Check if we found any valid amino acids
    if not composition:
        raise ProteinMetricsError(
            "Sequence contains no valid amino acid residues from the standard 20"
        )
    
    return composition


def molecular_weight(seq: str) -> float:
    """
    Calculate the molecular weight of a protein sequence.
    
    Computes the approximate molecular weight in Daltons (Da) based on the
    sum of average isotopic masses of constituent amino acids. The calculation
    accounts for peptide bond formation (loss of water molecules between residues).
    
    Formula: MW = sum(amino_acid_weights) + 18.015 (for terminal H and OH)
    
    Note: This is an approximation. Actual molecular weight may vary due to
    post-translational modifications, disulfide bonds, or prosthetic groups.
    Only the 20 standard amino acids contribute to the weight calculation.
    
    Args:
        seq (str): Protein sequence (whitespace is ignored)
        
    Returns:
        float: Molecular weight in Daltons (Da)
        
    Raises:
        ProteinMetricsError: If sequence is empty or contains no valid amino acids
                            from the standard 20
        
    Example:
        >>> molecular_weight('AAA')  # 3 alanines
        285.28
        >>> molecular_weight('GGG')  # 3 glycines
        243.22
    """
    cleaned = _clean_sequence(seq)
    
    if not cleaned:
        raise ProteinMetricsError("Sequence is empty or contains no valid characters")
    
    # Calculate sum of amino acid weights
    total_weight = 0.0
    valid_residues = 0
    
    for aa in cleaned:
        if aa in AMINO_ACID_WEIGHTS:
            total_weight += AMINO_ACID_WEIGHTS[aa]
            valid_residues += 1
    
    # Check if we found any valid amino acids
    if valid_residues == 0:
        raise ProteinMetricsError(
            "Sequence contains no valid amino acid residues from the standard 20"
        )
    
    # Add mass of terminal H (1.008) and OH (17.007) groups
    # These are lost during peptide bond formation for internal residues,
    # but remain at the N- and C-termini
    terminal_mass = 18.015  # H2O = 1.008 + 17.007
    
    return total_weight + terminal_mass


def get_summary_statistics(seq: str) -> Dict[str, Any]:
    """
    Get comprehensive summary statistics for a protein sequence.
    
    Convenience function that calculates all metrics at once.
    
    Args:
        seq (str): Protein sequence (whitespace is ignored)
        
    Returns:
        dict: Dictionary containing:
            - 'length': Sequence length (number of residues)
            - 'molecular_weight': Molecular weight in Daltons
            - 'composition': Amino acid composition dictionary
            - 'unique_amino_acids': Number of different amino acids present
            
    Raises:
        ProteinMetricsError: If sequence is empty or invalid
        
    Example:
        >>> stats = get_summary_statistics('MKLVVVGA')
        >>> stats['length']
        8
        >>> stats['molecular_weight']
        758.99
    """
    return {
        'length': length(seq),
        'molecular_weight': molecular_weight(seq),
        'composition': amino_acid_composition(seq),
        'unique_amino_acids': len(amino_acid_composition(seq))
    }


def get_amino_acid_percentage(seq: str) -> Dict[str, float]:
    """
    Calculate the percentage of each amino acid in the sequence.
    
    Helper function that provides composition as percentages rather than counts.
    Useful for comparing sequences of different lengths.
    
    Args:
        seq (str): Protein sequence (whitespace is ignored)
        
    Returns:
        dict: Dictionary mapping amino acid codes to their percentage (0.0-100.0)
              Only amino acids present in the sequence are included.
        
    Raises:
        ProteinMetricsError: If sequence is empty or invalid
        
    Example:
        >>> get_amino_acid_percentage('AAAGGG')
        {'A': 50.0, 'G': 50.0}
        >>> get_amino_acid_percentage('AAA')
        {'A': 100.0}
    """
    composition = amino_acid_composition(seq)
    seq_length = length(seq)
    
    # Calculate percentages for each amino acid
    percentages = {
        aa: (count / seq_length) * 100.0
        for aa, count in composition.items()
    }
    
    return percentages


if __name__ == "__main__":
    # Test/demonstration code
    print("Protein Metrics Module - Test Suite\n")
    print("=" * 80)
    
    test_sequences = [
        ("Short peptide", "MKLVVVGA"),
        ("Alanine repeat", "AAAAAAAAAA"),
        ("Mixed composition", "ACDEFGHIKLMNPQRSTVWY"),  # All 20 amino acids
        ("With whitespace", "MKLT VVGA\nGGKS"),
        ("Single residue", "M"),
        ("Glycine-rich", "GGGGGGGGGG"),
        ("Tryptophan (heaviest)", "WWWWW"),
    ]
    
    for name, seq in test_sequences:
        print(f"\nTest: {name}")
        print(f"Sequence: '{seq}'")
        
        try:
            stats = get_summary_statistics(seq)
            
            print(f"Length: {stats['length']} residues")
            print(f"Molecular Weight: {stats['molecular_weight']:.2f} Da")
            print(f"Unique amino acids: {stats['unique_amino_acids']}")
            
            # Show composition for shorter sequences
            if stats['length'] <= 20:
                comp_str = ', '.join(f"{aa}={count}" for aa, count in 
                                    sorted(stats['composition'].items()))
                print(f"Composition: {comp_str}")
            else:
                print(f"Composition: {len(stats['composition'])} different amino acids")
                
        except ProteinMetricsError as e:
            print(f"Error: {e}")
    
    print("\n" + "=" * 80)
    
    # Test edge cases
    print("\nEdge Case Tests:")
    print("-" * 80)
    
    edge_cases = [
        ("Empty string", ""),
        ("Only whitespace", "   \n\t  "),
        ("Non-standard letters (X,Y,Z)", "XYZ123"),
    ]
    
    for name, seq in edge_cases:
        print(f"\n{name}: '{seq}'")
        try:
            result = length(seq)
            print(f"  Length: {result}")
        except ProteinMetricsError as e:
            print(f"  Error (expected): {e}")
    
    print("\n" + "=" * 80)
    print("\nTo use this module, import the functions:")
    print("  from protein_metrics import amino_acid_composition, molecular_weight, length")
