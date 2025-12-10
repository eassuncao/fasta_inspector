#!/usr/bin/env python3
# AI Assistance: GitHub Copilot used for scaffolding and function orchestration.

"""
FASTA Inspector - Bioinformatics Sequence Analysis Tool

This program analyzes multi-FASTA files containing DNA and protein sequences.
Sequences classified as DNA may actually be DNA or RNA - the tool treats U (uracil)
as equivalent to T (thymine) for classification and metrics calculations.

For each sequence, the program automatically classifies the type and computes:
  - DNA/RNA sequences: GC content, AT/GC ratio, N count, base composition
    (Note: RNA bases with U are mapped to the T bucket for base counts)
  - Protein sequences: molecular weight, amino acid composition, length

Output includes:
  - Per-sequence reports with type, header, and relevant metrics
  - Final summary showing total sequences analyzed and counts/percentages for
    DNA/RNA vs protein sequences

The tool demonstrates a modular architecture with separate modules for:
  - FASTA file parsing (fasta_reader)
  - Sequence type classification (type_classifier)
  - DNA metrics calculation (dna_metrics)
  - Protein metrics calculation (protein_metrics)

Usage:
    python main.py [fasta_file]
    
If no file is specified, defaults to 'test_mixed.fasta' in the current directory.
"""

import sys
import os
import time
from typing import Tuple

# Import project modules
from fasta_reader import read_fasta, FastaFormatError
from type_classifier import detect_sequence_type, SequenceClassificationError
from dna_metrics import (
    gc_content, 
    n_count, 
    at_gc_ratio, 
    length as dna_length,
    get_base_counts
)
from protein_metrics import (
    amino_acid_composition,
    molecular_weight,
    length as protein_length,
    ProteinMetricsError
)
from utils import format_large_number, truncate_string


def print_header():
    """Print the program header."""
    print("=" * 80)
    print("FASTA INSPECTOR - Sequence Analysis Tool".center(80))
    print("=" * 80)
    print()


def print_separator():
    """Print a section separator."""
    print("-" * 80)


def analyze_dna_sequence(sequence: str, seq_id: int, header: str) -> None:
    """
    Analyze and report metrics for a DNA or RNA sequence.
    
    RNA sequences containing U (uracil) are supported - U is treated as
    equivalent to T for all metrics calculations.
    
    Args:
        sequence (str): The DNA or RNA sequence to analyze
        seq_id (int): Sequence number in the file
        header (str): FASTA header line
    """
    print(f"\nSequence #{seq_id}: DNA/RNA")
    print_separator()
    print(f"Header: {truncate_string(header, 70)}")
    print()
    
    try:
        # Calculate DNA metrics
        seq_length = dna_length(sequence)
        gc_pct = gc_content(sequence)
        n_bases = n_count(sequence)
        at_gc = at_gc_ratio(sequence)
        base_comp = get_base_counts(sequence)
        
        # Print metrics
        print(f"Length:              {format_large_number(seq_length)} bp")
        print(f"GC Content:          {gc_pct:.2f}%")
        print(f"AT/GC Ratio:         {at_gc:.4f}" if at_gc != float('inf') else "AT/GC Ratio:         inf (no GC bases)")
        print(f"Ambiguous Bases (N): {n_bases}")
        
        # Base composition
        # Note: T includes U for RNA sequences. The dna_metrics module maps U→T.
        # "DNA/RNA" label reflects that sequences with U are allowed and handled consistently.
        print(f"\nBase Composition:")
        print(f"  A (Adenine):       {base_comp['A']:>6}  ({base_comp['A']/seq_length*100:5.2f}%)")
        print(f"  T (Thymine/U):     {base_comp['T']:>6}  ({base_comp['T']/seq_length*100:5.2f}%)")
        print(f"  C (Cytosine):      {base_comp['C']:>6}  ({base_comp['C']/seq_length*100:5.2f}%)")
        print(f"  G (Guanine):       {base_comp['G']:>6}  ({base_comp['G']/seq_length*100:5.2f}%)")
        
        if base_comp['N'] > 0:
            print(f"  N (Ambiguous):     {base_comp['N']:>6}  ({base_comp['N']/seq_length*100:5.2f}%)")
        
        if base_comp['other'] > 0:
            print(f"  Other codes:       {base_comp['other']:>6}  ({base_comp['other']/seq_length*100:5.2f}%)")
        
        print()
        
    except Exception as e:
        print(f"Error analyzing DNA sequence: {e}")
        print()


def analyze_protein_sequence(sequence: str, seq_id: int, header: str) -> None:
    """
    Analyze and report metrics for a protein sequence.
    
    Prints length (in residues), molecular weight (in Da and kDa), number of
    unique amino acids present, and the top 10 most abundant amino acids with
    their counts and percentages.
    
    Args:
        sequence (str): The protein sequence to analyze
        seq_id (int): Sequence number in the file
        header (str): FASTA header line
    """
    print(f"\nSequence #{seq_id}: PROTEIN")
    print_separator()
    print(f"Header: {truncate_string(header, 70)}")
    print()
    
    try:
        # Calculate protein metrics
        seq_length = protein_length(sequence)
        mol_weight = molecular_weight(sequence)
        aa_comp = amino_acid_composition(sequence)
        
        # Print metrics
        print(f"Length:              {format_large_number(seq_length)} residues")
        print(f"Molecular Weight:    {format_large_number(mol_weight, 2)} Da ({mol_weight/1000:.2f} kDa)")
        print(f"Unique Amino Acids:  {len(aa_comp)}/20")
        
        # Amino acid composition - show most abundant
        print(f"\nTop 10 Most Abundant Amino Acids:")
        sorted_aa = sorted(aa_comp.items(), key=lambda x: x[1], reverse=True)[:10]
        
        for aa, count in sorted_aa:
            percentage = (count / seq_length) * 100
            print(f"  {aa}: {count:>4}  ({percentage:5.2f}%)")
        
        print()
        
    except ProteinMetricsError as e:
        print(f"Error analyzing protein sequence: {e}")
        print()
    except Exception as e:
        print(f"Unexpected error: {e}")
        print()


def analyze_fasta_file(filepath: str) -> Tuple[int, int, int, int, int]:
    """
    Analyze all sequences in a FASTA file.
    
    Iterates over all sequences using read_fasta, classifies each sequence
    with detect_sequence_type, and routes to the appropriate analysis function
    (DNA/RNA or protein). Catches and reports FASTA format errors and
    classification errors while continuing with remaining sequences.
    
    Args:
        filepath (str): Path to the FASTA file
        
    Returns:
        tuple: (total_sequences, dna_count, protein_count, total_bases, total_residues)
    """
    total_sequences = 0
    dna_count = 0
    protein_count = 0
    total_bases = 0
    total_residues = 0
    
    try:
        for seq_id, (header, sequence) in enumerate(read_fasta(filepath), 1):
            total_sequences += 1
            
            try:
                # Classify the sequence
                seq_type = detect_sequence_type(sequence)
                
                # Route to appropriate analysis function
                if seq_type == 'DNA':
                    analyze_dna_sequence(sequence, seq_id, header)
                    dna_count += 1
                    total_bases += dna_length(sequence)
                elif seq_type == 'PROTEIN':
                    analyze_protein_sequence(sequence, seq_id, header)
                    protein_count += 1
                    total_residues += protein_length(sequence)
                else:
                    print(f"\nSequence #{seq_id}: UNKNOWN TYPE")
                    print_separator()
                    print(f"Header: {truncate_string(header, 70)}")
                    print(f"Could not classify sequence type: {seq_type}")
                    print()
                    
            except SequenceClassificationError as e:
                print(f"\nSequence #{seq_id}: CLASSIFICATION ERROR")
                print_separator()
                print(f"Header: {truncate_string(header, 70)}")
                print(f"Error: {e}")
                print()
                
    except FastaFormatError as e:
        print(f"\nFASTA Format Error: {e}")
        return total_sequences, dna_count, protein_count, total_bases, total_residues
    except FileNotFoundError as e:
        print(f"\nFile Error: {e}")
        return total_sequences, dna_count, protein_count, total_bases, total_residues
    except Exception as e:
        print(f"\nUnexpected Error: {e}")
        return total_sequences, dna_count, protein_count, total_bases, total_residues
    
    return total_sequences, dna_count, protein_count, total_bases, total_residues


def print_summary(total: int, dna: int, protein: int):
    """
    Print analysis summary.
    
    Displays total number of sequences analyzed, count and percentage of
    DNA/RNA sequences, count and percentage of protein sequences, and if
    applicable, the number and percentage of unclassified sequences.
    
    Args:
        total (int): Total number of sequences
        dna (int): Number of DNA/RNA sequences
        protein (int): Number of protein sequences
    """
    print("=" * 80)
    print("ANALYSIS SUMMARY".center(80))
    print("=" * 80)
    print()
    print(f"Total Sequences:     {total}")
    print(f"DNA/RNA Sequences:   {dna} ({dna/total*100:.1f}%)" if total > 0 else "DNA/RNA Sequences:   0")
    print(f"Protein Sequences:   {protein} ({protein/total*100:.1f}%)" if total > 0 else "Protein Sequences:   0")
    
    unclassified = total - dna - protein
    if unclassified > 0:
        print(f"Unclassified:        {unclassified} ({unclassified/total*100:.1f}%)")
    
    print()
    print("=" * 80)


def main():
    """
    Main entry point for the FASTA Inspector program.
    
    Handles command-line argument for the FASTA file path (or defaults to
    'test_mixed.fasta' if no argument is provided). Prints basic file
    information (name and size), calls analyze_fasta_file to process all
    sequences, and prints a summary of the results.
    """
    start_time = time.perf_counter()
    
    print_header()
    
    # Determine input file
    if len(sys.argv) > 1:
        fasta_file = sys.argv[1]
    else:
        fasta_file = "test_mixed.fasta"
    
    # Check if file exists
    if not os.path.isfile(fasta_file):
        print(f"Error: File '{fasta_file}' not found.")
        print()
        print("Usage: python main.py [fasta_file]")
        print("       If no file is specified, the program defaults to 'test_mixed.fasta'")
        print()
        sys.exit(1)
    
    # Display file info
    file_size = os.path.getsize(fasta_file)
    print(f"Analyzing file: {fasta_file}")
    print(f"File size:      {file_size:,} bytes ({file_size/1024:.2f} KB)")
    print()
    
    # Analyze the file
    total, dna, protein, total_bases, total_residues = analyze_fasta_file(fasta_file)
    
    # Print benchmarking results
    end_time = time.perf_counter()
    elapsed = end_time - start_time
    
    print(f"\nProcessed {total} sequences in {elapsed:.2f} seconds")
    rate = total / elapsed if elapsed > 0 else 0
    print(f"≈ {rate:,.0f} sequences/second")
    
    print(f"\nTotal DNA bases processed: {total_bases:,}")
    print(f"Total protein residues processed: {total_residues:,}")
    print()
    
    # Print summary
    if total > 0:
        print_summary(total, dna, protein)
    else:
        print("No sequences were successfully analyzed.")
        print()


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\nAnalysis interrupted by user.")
        print()
        sys.exit(0)
    except Exception as e:
        print(f"\n\nFatal error: {e}")
        print()
        sys.exit(1)
