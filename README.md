A complete Python bioinformatics pipeline that combines DNA analysis,
transcription, translation, mutation detection and generates a
professional HTML report with dashboard charts.
Built as Phase 1 - Project 5 of my bioinformatics-to-tech portfolio.

## What it does
Runs 4 modules automatically in one command:
- Module 1 - DNA sequence analysis (GC content, base count, complement)
- Module 2 - DNA to RNA transcription and protein translation
- Module 3 - Mutation detection (SNPs, transitions, transversions)
- Module 4 - HTML report with dashboard charts

## Features
- Full automated pipeline with progress bar
- Command line arguments support
- Organized output folder with all files
- Professional HTML report with statistics dashboard
- 6 panel visualization chart
- Processes any DNA sequence in seconds

## Sample Output
Sequence length    : 33 bases
GC content         : 45.45%
Codons translated  : 6
STOP codon found   : Yes
Mutations detected : 3
Mutation rate      : 9.09%
Similarity         : 90.91%

## Usage

Interactive mode:
python pipeline.py

Command line mode:
python pipeline.py --sequence ATGCTTGAATTTGCCTAA --reference ATGCAAGAATTTGCCTAA

## Tech Stack
- Python 3.13
- Matplotlib
- VS Code
- Git + GitHub

## Project Roadmap
- Phase 1 Project 1 - DNA Sequence Analyzer - Done
- Phase 1 Project 2 - DNA to RNA Transcriber - Done
- Phase 1 Project 3 - FASTA File Parser - Done
- Phase 1 Project 4 - Mutation Detector - Done
- Phase 1 Project 5 - Mini Pipeline - Done
- Phase 2 - Intermediate Bioinformatics - Coming soon

## Author
Padma Shree Jena
Bioinformatics + Tech Enthusiast | Python | R | Bash
GitHub: https://github.com/Paddu2006

## License
MIT License
