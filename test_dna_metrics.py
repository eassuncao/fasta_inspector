#!/usr/bin/env python3
"""Quick test of dna_metrics with actual FASTA sequences."""

from fasta_reader import read_fasta
from type_classifier import detect_sequence_type
from dna_metrics import gc_content, n_count, at_gc_ratio, length, get_summary_statistics

# Read test file
print("Testing DNA Metrics with test_mixed.fasta\n")
print("=" * 80)

for i, (header, sequence) in enumerate(read_fasta('test_mixed.fasta'), 1):
    seq_type = detect_sequence_type(sequence)
    
    # Only analyze DNA sequences
    if seq_type == 'DNA':
        print(f"\nSequence {i}: DNA")
        print(f"Header: {header[:60]}...")
        
        stats = get_summary_statistics(sequence)
        
        print(f"Length: {stats['length']} bp")
        print(f"GC Content: {stats['gc_content']:.2f}%")
        print(f"N Count: {stats['n_count']}")
        
        ratio = stats['at_gc_ratio']
        if ratio == float('inf'):
            print(f"AT/GC Ratio: inf")
        else:
            print(f"AT/GC Ratio: {ratio:.4f}")
        
        print(f"Composition: A={stats['base_counts']['A']}, "
              f"T={stats['base_counts']['T']}, "
              f"C={stats['base_counts']['C']}, "
              f"G={stats['base_counts']['G']}")

print("\n" + "=" * 80)
