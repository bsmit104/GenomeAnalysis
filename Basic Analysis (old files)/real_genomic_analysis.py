# Brayden Smith
# jbrayden35@gmail.com
# script for analyzing k mer frequencies in real genome

import time
import random
from collections import defaultdict
import multiprocessing as mp
import matplotlib.pyplot as plt
import argparse
import os

def read_fasta(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return sequence

def kmer_frequency_basic(genome, k):
    start_time = time.time()
    kmer_counts = defaultdict(int)
    for i in range(len(genome) - k + 1):
        kmer = genome[i:i + k]
        kmer_counts[kmer] += 1
    end_time = time.time()
    return dict(kmer_counts), end_time - start_time

def count_kmers_chunk(chunk, k):
    kmer_counts = defaultdict(int)
    for i in range(len(chunk) - k + 1):
        kmer = chunk[i:i + k]
        kmer_counts[kmer] += 1
    return dict(kmer_counts)

def combine_counts(results):
    combined = defaultdict(int)
    for d in results:
        for kmer, count in d.items():
            combined[kmer] += count
    return dict(combined)

def kmer_frequency_parallel(genome, k):
    start_time = time.time()
    num_cores = mp.cpu_count()
    chunk_size = len(genome) // num_cores
    chunks = [genome[i:i + chunk_size + k - 1] for i in range(0, len(genome), chunk_size)]
    if len(chunks[-1]) < k:
        chunks[-1] = genome[-chunk_size - k + 1:]
    
    with mp.Pool(num_cores) as pool:
        results = pool.starmap(count_kmers_chunk, [(chunk, k) for chunk in chunks])
    
    counts = combine_counts(results)
    end_time = time.time()
    return counts, end_time - start_time

def main():
    parser = argparse.ArgumentParser(description="Analyze k-mer frequencies in a genomic FASTA file.")
    parser.add_argument("filename", help="Name of the FASTA file in 'real_genomic_data' folder (ex. sequence.fasta)")
    args = parser.parse_args()

    # Construct file path
    data_folder = "real_genomic_data"
    file_path = os.path.join(data_folder, args.filename)

    if not os.path.exists(file_path):
        print(f"Error: File '{file_path}' not found. Make sure it's in the '{data_folder}' folder.")
        return

    genome = read_fasta(file_path)
    print(f"Loaded '{args.filename}' with {len(genome)} bases")
    print(f"First 50 bases: {genome[:50]}")
    
    k = 5  # 5-mers
    
    print("Running comparisons...")
    counts_b, time_b = kmer_frequency_basic(genome, k)
    counts_p, time_p = kmer_frequency_parallel(genome, k)
    
    print(f"Basic: {time_b:.4f} sec, Unique k-mers: {len(counts_b)}")
    print(f"Parallel: {time_p:.4f} sec, Unique k-mers: {len(counts_p)}")
    speedup = (time_b - time_p) / time_b * 100
    print(f"Speedup: {speedup:.2f}%")


    output_file = f"kmer_results_{args.filename.split('.')[0]}.txt"
    with open(output_file, "w") as f:
        f.write(f"File: {args.filename}\n")
        f.write(f"Basic: {time_b:.4f} sec\nParallel: {time_p:.4f} sec\nSpeedup: {speedup:.2f}%\n")
        f.write(f"Top 10 {k}-mers:\n")
        for kmer, count in sorted(counts_p.items(), key=lambda x: x[1], reverse=True)[:10]:
            f.write(f"{kmer}: {count}\n")
    

    top_kmers = sorted(counts_p.items(), key=lambda x: x[1], reverse=True)[:10]
    kmers, freqs = zip(*top_kmers)
    plt.bar(kmers, freqs)
    plt.xlabel(f"{k}-mers")
    plt.ylabel("Frequency")
    plt.title(f"Top 10 {k}-mers in {args.filename.split('.')[0]}")
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f"kmer_plot_{args.filename.split('.')[0]}.png")
    print(f"Results saved to '{output_file}' and 'kmer_plot_{args.filename.split('.')[0]}.png'")

if __name__ == "__main__":
    main()