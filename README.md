# FASTA Inspector

A bioinformatics tool for analyzing FASTA sequence files, automatically classifying sequences as DNA/RNA or protein, and computing relevant metrics for each sequence type.

## Project Description

FASTA Inspector is a command-line utility designed to parse FASTA files, detect sequence types, and calculate comprehensive statistics. The tool handles mixed FASTA files containing both nucleotide and protein sequences, applying appropriate analyses based on automatic sequence classification.

## Dataset Summary

The project includes test datasets representing diverse biological sequences:

- **E. coli CDS DNA**: Coding sequence nucleotide data from *Escherichia coli*
- **E. coli Protein CDS**: Translated protein sequences from *E. coli* coding regions
- **HLA-B27 Nucleotide**: Human leukocyte antigen B27 gene sequences
- **HLA-B27 Protein**: HLA-B27 protein sequences
- **Human Complement Protein**: Complement system protein sequences from UniProt
- **test_mixed.fasta**: Combined testing file with 13 sequences (4 DNA, 9 protein) for validation

## Features

### Core Functionality
- **Automatic Sequence Type Detection**: Classifies sequences as DNA/RNA or protein using nucleotide composition thresholds
- **FASTA Parsing**: Robust parser with fallback header detection for non-standard formats
- **Mixed File Support**: Handles FASTA files containing both nucleotide and protein sequences

### DNA/RNA Metrics
- Sequence length
- GC content percentage
- AT/GC ratio
- Ambiguous nucleotide (N) count
- Base composition (A, T/U, G, C counts)

### Protein Metrics
- Sequence length
- Amino acid composition (all 20 standard amino acids)
- Molecular weight calculation
- Per-residue percentage distribution

### Output
- Summary statistics for each sequence
- Sequence type classification
- Comprehensive metrics formatted for readability
- Classification summary (total sequences, DNA count, protein count)

## How to Run

### Basic Usage
```bash
python main.py
```

The program will prompt you to enter the path to a FASTA file. Alternatively, run with a file argument:

```bash
python main.py path/to/your/file.fasta
```

### Example Output
```
Enter the path to the FASTA file: test_mixed.fasta

Processing FASTA file: test_mixed.fasta
...
[Sequence analysis results]
...

=== Classification Summary ===
Total sequences: 13
DNA sequences: 4 (30.8%)
Protein sequences: 9 (69.2%)
```

## Requirements

- **Python**: 3.10 or higher
- **External Libraries**: None required (uses Python standard library only)

### Development Dependencies
For running the test suite:
```bash
pip install pytest>=7.0.0
```

### Running Tests
```bash
pytest -v
```

The project includes 188 comprehensive tests covering all modules.

## Development Notes

### AI Assistance
This project was developed with AI assistance (GitHub Copilot) for:
- Code implementation and documentation
- Test suite creation (pytest-based)
- Code refactoring and optimization
- Documentation generation

### Version Control
The project uses Git for version control. Key development milestones:
- Initial module implementation
- Comprehensive pytest test suite (188 tests)
- Documentation and code cleanup
- README finalization

## Limitations and Next Steps

### Current Limitations
- **Classification Threshold**: Uses a fixed 90% nucleotide composition threshold; ambiguous sequences (85-95%) may be misclassified
- **Standard Alphabets Only**: Limited support for non-standard IUPAC ambiguity codes beyond common bases
- **Single File Processing**: Processes one file at a time; no batch processing capability
- **Output Format**: Console output only; no export to CSV/JSON
- **Memory Usage**: Loads entire file into memory; large files (>1GB) may cause issues

### Planned Enhancements
- **Configurable Thresholds**: Allow user-specified classification thresholds
- **Batch Processing**: Process multiple FASTA files in a directory
- **Export Options**: Add CSV, JSON, and TSV output formats
- **Visualization**: Generate plots for composition analysis and sequence statistics
- **Advanced Metrics**: 
  - Codon usage analysis for DNA
  - Hydrophobicity profiles for proteins
  - ORF detection and translation
- **Streaming Parser**: Handle large files without loading entirely into memory
- **Web Interface**: Optional GUI for non-command-line users

## Project Structure

```
fasta_inspector/
├── main.py                    # Entry point and orchestration
├── fasta_reader.py            # FASTA file parsing
├── type_classifier.py         # Sequence type detection
├── dna_metrics.py             # DNA/RNA analysis
├── protein_metrics.py         # Protein analysis
├── utils.py                   # Shared utilities
├── test_*.py                  # Pytest test suite (6 modules)
├── test_mixed.fasta           # Test dataset
├── requirements.txt           # Python dependencies
└── README.md                  # This file
```

## License

This project was created for academic purposes (COMP625 - Algorithms for Bioinformatics, Athabasca University).

## Author

Erich Assunção
Athabasca University  
December 2025
