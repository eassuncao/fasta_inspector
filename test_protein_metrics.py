#!/usr/bin/env python3
"""Quick test of protein_metrics with actual FASTA sequences."""

from fasta_reader import read_fasta
from type_classifier import detect_sequence_type
from protein_metrics import (
    amino_acid_composition, 
    molecular_weight, 
    length, 
    get_summary_statistics
)

# Read test file
print("Testing Protein Metrics with test_mixed.fasta\n")
print("=" * 80)

for i, (header, sequence) in enumerate(read_fasta('test_mixed.fasta'), 1):
    seq_type = detect_sequence_type(sequence)
    
    # Only analyze protein sequences
    if seq_type == 'PROTEIN':
        print(f"\nSequence {i}: PROTEIN")
        print(f"Header: {header[:60]}...")
        
        stats = get_summary_statistics(sequence)
        
        print(f"Length: {stats['length']} residues")
        print(f"Molecular Weight: {stats['molecular_weight']:.2f} Da ({stats['molecular_weight']/1000:.2f} kDa)")
        print(f"Unique amino acids: {stats['unique_amino_acids']}/20")
        
        # Show top 5 most abundant amino acids
        sorted_comp = sorted(stats['composition'].items(), 
                            key=lambda x: x[1], reverse=True)
        top_5 = sorted_comp[:5]
        top_str = ', '.join(f"{aa}={count}" for aa, count in top_5)
        print(f"Most abundant: {top_str}")

print("\n" + "=" * 80)
