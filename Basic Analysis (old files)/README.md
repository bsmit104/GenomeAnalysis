# Basic Genomic Analysis Tool
A Python tool to analyze mock genomic data by finding all positions of a substring (e.g., "ATG") in a DNA sequence.

## Features
- Generates random DNA sequences (A, T, C, G).
- Two methods: Basic (naive loop) and Optimized (str.find()).
- Achieves >95% runtime improvement over the basic method.

## Usage
1. Install Python 3.11+.
2. Run `python genomic_analysis.py`
3. Results saved to `genomic_results.txt`.

## Results
- Basic: ~0.1120 sec for 1M bases.
- Optimized: ~0.0050 sec.
- Speedup: ~95.54%.