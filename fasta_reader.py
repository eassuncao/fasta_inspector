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
                
                # Header line (starts with '>')
                if line.startswith('>'):
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
                    current_header = line[1:].strip()  # Remove '>' and strip whitespace
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
    
    Args:
        filepath (str): Path to the file to validate
        
    Returns:
        bool: True if the file is valid, False otherwise
        
    Example:
        >>> if validate_fasta_file('sequences.fasta'):
        ...     print("File is valid")
    """
    try:
        # Try to read all sequences
        count = sum(1 for _ in read_fasta(filepath))
        return count > 0
    except (FileNotFoundError, FastaFormatError, IOError):
        return False


def count_sequences(filepath: str) -> int:
    """
    Count the number of sequences in a FASTA file.
    
    Args:
        filepath (str): Path to the FASTA file
        
    Returns:
        int: Number of sequences in the file
        
    Raises:
        FileNotFoundError: If the file does not exist
        FastaFormatError: If the file format is invalid
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
