# AI Assistance: GitHub Copilot helpful in generating this module.

"""
FASTA File Reader Module

This module provides functionality for reading and parsing FASTA format files.
It supports multi-FASTA files and handles multiline sequences automatically.

Functions:
    read_fasta(filepath): Generator that yields (header, sequence) tuples
"""

import os
from typing import Generator, Tuple


class FastaFormatError(Exception):
    """Exception raised for errors in FASTA file format."""
    pass


def is_header_line(line: str) -> bool:
    """
    Determine if a line should be treated as a FASTA header.
    
    Returns True if:
    1. Line starts with '>' (standard FASTA header), OR
    2. Line appears to be a header-like metadata line (fallback detection)
    
    Fallback rule: A line is a header if it doesn't look like sequence data.
    Sequence data is typically:
    - Longer stretches of characters (>60 chars typically)
    - Contains mostly biological sequence characters
    - Lacks formatting characters like |, :, [, ], spaces in structured patterns
    
    Args:
        line (str): The line to check (should already be stripped)
        
    Returns:
        bool: True if this line should be treated as a header
        
    Example:
        >>> is_header_line('>seq1 description')
        True
        >>> is_header_line('ATCGATCGATCG')
        False
        >>> is_header_line('sp|P12345|PROTEIN_NAME')
        True
    """
    # Validate input
    if not line or not isinstance(line, str):
        return False
    
    # Standard FASTA header starts with '>'
    if line.startswith('>'):
        return True
    
    # Fallback detection for non-standard headers
    # Check for common header formatting patterns
    header_indicators = ['|', ':', '[', ']', 'sp|', 'tr|', 'gb|', 'ref|', 
                        'pdb|', 'GenBank', 'UniProt', 'Accession', 'lcl|']
    if any(indicator in line for indicator in header_indicators):
        return True
    
    # Count alphabetic characters
    alpha_chars = [c for c in line.upper() if c.isalpha()]
    if not alpha_chars:
        # Line has no alphabetic characters - likely a header (e.g., numeric IDs)
        return True
    
    # Define sequence character sets
    # DNA/RNA bases including ambiguity codes
    dna_chars = set('ACGTUNRYSWKMBDHV')
    # Standard 20 amino acids
    protein_chars = set('ACDEFGHIKLMNPQRSTVWY')
    # Combined sequence alphabet
    sequence_chars = dna_chars | protein_chars
    
    # Count how many characters are valid sequence characters
    seq_char_count = sum(1 for c in alpha_chars if c in sequence_chars)
    
    # Calculate proportion of sequence characters
    if len(alpha_chars) > 0:
        seq_proportion = seq_char_count / len(alpha_chars)
        
        # If >85% of alphabetic chars are sequence chars, it's likely sequence data
        if seq_proportion > 0.85:
            return False
    
    # Otherwise, treat as a header (likely metadata, accession ID, etc.)
    return True


def read_fasta(filepath: str) -> Generator[Tuple[str, str], None, None]:
    """
    Read and parse a FASTA format file.
    
    This function reads a FASTA file and yields each sequence entry as a 
    (header, sequence) tuple. It supports multi-FASTA files with multiple 
    sequences and automatically handles multiline sequences by joining them.
    
    Args:
        filepath (str): Path to the FASTA file to read
        
    Yields:
        tuple: A (header, sequence) tuple where:
            - header (str): The header line without the '>' character
            - sequence (str): The complete sequence (multiline sequences are joined)
            
    Raises:
        FileNotFoundError: If the specified file does not exist
        FastaFormatError: If the file format is invalid (e.g., no header, 
                         sequence before header, empty sequences)
        
    Example:
        >>> for header, sequence in read_fasta('sequences.fasta'):
        ...     print(f"Header: {header}")
        ...     print(f"Sequence length: {len(sequence)}")
    """
    # Check if file exists
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"FASTA file not found: {filepath}")
    
    # Check if file is readable
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"Path is not a file: {filepath}")
    
    current_header = None
    current_sequence = []
    line_number = 0
    has_content = False
    
    try:
        with open(filepath, 'r', encoding='utf-8') as file:
            for line in file:
                line_number += 1
                line = line.strip()
                
                # Skip blank lines
                if not line:
                    continue
                
                has_content = True
                
                # Check if this line is a header (either standard '>' or fallback detection)
                if is_header_line(line):
                    # If we have a previous sequence, yield it
                    if current_header is not None:
                        # Join all sequence lines and validate
                        sequence = ''.join(current_sequence)
                        if not sequence:
                            raise FastaFormatError(
                                f"Empty sequence for header '{current_header}' "
                                f"at line {line_number}"
                            )
                        yield (current_header, sequence)
                    
                    # Start new sequence
                    # Remove '>' if present, otherwise use the line as-is
                    if line.startswith('>'):
                        current_header = line[1:].strip()
                    else:
                        current_header = line.strip()
                    
                    if not current_header:
                        raise FastaFormatError(
                            f"Empty header at line {line_number}"
                        )
                    current_sequence = []
                
                # Sequence line
                else:
                    # Make sure we have a header before sequence data
                    if current_header is None:
                        raise FastaFormatError(
                            f"Sequence data found before header at line {line_number}: '{line}'"
                        )
                    # Add sequence line (already stripped of whitespace)
                    current_sequence.append(line)
            
            # Yield the last sequence if present
            if current_header is not None:
                sequence = ''.join(current_sequence)
                if not sequence:
                    raise FastaFormatError(
                        f"Empty sequence for header '{current_header}' "
                        f"at end of file"
                    )
                yield (current_header, sequence)
            elif has_content:
                # File had content but no valid sequences
                raise FastaFormatError(
                    "File contains data but no valid FASTA sequences found"
                )
            else:
                # Empty file
                raise FastaFormatError("File is empty")
                
    except UnicodeDecodeError as e:
        raise FastaFormatError(
            f"File encoding error at line {line_number}: {str(e)}"
        )
    except IOError as e:
        raise IOError(f"Error reading file '{filepath}': {str(e)}")


def validate_fasta_file(filepath: str) -> bool:
    """
    Validate that a file is a properly formatted FASTA file.
    
    This function attempts to parse the file and returns True if at least
    one valid sequence is found. It catches all exceptions and returns False
    for any invalid files.
    
    Args:
        filepath (str): Path to the file to validate
        
    Returns:
        bool: True if the file is valid and contains at least one sequence,
              False otherwise
        
    Example:
        >>> if validate_fasta_file('sequences.fasta'):
        ...     print("File is valid")
        ...     # Proceed with processing
    """
    try:
        # Try to read all sequences - if this succeeds, file is valid
        count = sum(1 for _ in read_fasta(filepath))
        return count > 0
    except (FileNotFoundError, FastaFormatError, IOError):
        return False


def count_sequences(filepath: str) -> int:
    """
    Count the number of sequences in a FASTA file.
    
    This function parses the entire file to count sequences. For faster
    counting on large files (without validation), consider using
    utils.count_sequences_quick() instead.
    
    Args:
        filepath (str): Path to the FASTA file
        
    Returns:
        int: Number of sequences in the file
        
    Raises:
        FileNotFoundError: If the file does not exist
        FastaFormatError: If the file format is invalid
        
    Example:
        >>> count = count_sequences('sequences.fasta')
        >>> print(f"Found {count} sequences")
    """
    return sum(1 for _ in read_fasta(filepath))


if __name__ == "__main__":
    # Simple test/demonstration
    import sys
    
    # Default test file if no arguments provided
    if len(sys.argv) == 1:
        test_file = "test_mixed.fasta"  # Default file for VS Code Run button
        print(f"No file specified, using default: {test_file}\n")
    else:
        test_file = sys.argv[1]
    
    try:
        print(f"Reading FASTA file: {test_file}\n")
        for i, (header, sequence) in enumerate(read_fasta(test_file), 1):
            print(f"Sequence {i}:")
            print(f"  Header: {header}")
            print(f"  Length: {len(sequence)}")
            print(f"  Preview: {sequence[:50]}{'...' if len(sequence) > 50 else ''}")
            print()
    except (FileNotFoundError, FastaFormatError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
