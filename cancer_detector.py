import time
from collections import defaultdict
import multiprocessing as mp
import matplotlib.pyplot as plt
import argparse
import os
import pandas as pd
import numpy as np
from Bio import SeqIO, Seq

def read_fasta(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    sequence = ''.join(line.strip() for line in lines if not line.startswith('>'))
    return sequence

def load_gene_annotations(annotation_file, genome_type="gff3"):
    genes = {}
    df = pd.read_csv(annotation_file, sep='\t', comment='#', header=None, dtype=str)
    print(f"GFF3 seqids: {sorted(set(df[0].unique()))}")
    for _, row in df[df[2] == 'gene'].iterrows():
        gene_info = dict(item.split('=') for item in row[8].split(';'))
        gene_name = gene_info.get('Name', 'unknown')
        start, end = int(row[3]) - 1, int(row[4])
        genes[gene_name] = (start, end)
    return genes

def load_cancer_genes(cancer_file):
    df = pd.read_csv(cancer_file)
    return set(df['Gene Symbol'].dropna())

def kmer_frequency_basic(genome, k):
    start_time = time.time()
    kmer_counts = defaultdict(int)
    for i in range(len(genome) - k + 1):
        kmer = genome[i:i + k]
        if 'N' not in kmer:
            kmer_counts[kmer] += 1
    end_time = time.time()
    return dict(kmer_counts), end_time - start_time

def count_kmers_chunk(chunk, k):
    kmer_counts = defaultdict(int)
    for i in range(len(chunk) - k + 1):
        kmer = chunk[i:i + k]
        if 'N' not in kmer:
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

def find_cancer_kmers(genome, k, genes, cancer_genes):
    cancer_kmer_counts = defaultdict(int)
    cancer_gene_hits = defaultdict(list)
    gene_kmer_counts = defaultdict(lambda: defaultdict(int))  # k-mer counts per gene
    cancer_regions = set()
    for gene, (start, end) in genes.items():
        if gene in cancer_genes:
            gene_seq = genome[start:end]
            if len(gene_seq) >= k:
                for i in range(len(gene_seq) - k + 1):
                    kmer = gene_seq[i:i + k]
                    if 'N' not in kmer:
                        cancer_kmer_counts[kmer] += 1
                        cancer_gene_hits[kmer].append(gene)
                        gene_kmer_counts[gene][kmer] += 1  # Track per gene
                cancer_regions.update(range(start, end))
    total_cancer_bases = len(cancer_regions)
    print(f"Total bases in cancer genes: {total_cancer_bases}")
    return dict(cancer_kmer_counts), dict(cancer_gene_hits), cancer_regions, gene_kmer_counts

def main():
    parser = argparse.ArgumentParser(description="Analyze k-mers and cancer genes in a FASTA file.")
    parser.add_argument("filename", help="FASTA file in 'real_genomic_data' (e.g., human_chrom1.fa)")
    args = parser.parse_args()

    data_folder = "real_genomic_data"
    file_path = os.path.join(data_folder, args.filename)
    
    if not os.path.exists(file_path):
        print(f"Error: File '{file_path}' not found.")
        return

    genome = read_fasta(file_path)
    print(f"Loaded '{args.filename}' with {len(genome)} bases")
    print(f"First 50 bases: {genome[:50]}")

    annotation_file = os.path.join(data_folder, "human_chrom1.gff3")
    genes = load_gene_annotations(annotation_file, genome_type="gff3")
    cancer_genes = load_cancer_genes("cosmic_cancer_genes.csv")

    k = 5
    print("Running k-mer analysis...")
    counts_b, time_b = kmer_frequency_basic(genome, k)
    counts_p, time_p = kmer_frequency_parallel(genome, k)
    
    print(f"Basic: {time_b:.4f} sec, Unique k-mers: {len(counts_b)}")
    print(f"Parallel: {time_p:.4f} sec, Unique k-mers: {len(counts_p)}")
    speedup = (time_b - time_p) / time_b * 100
    print(f"Speedup: {speedup:.2f}%")

    # Find cancer k-mers with gene associations
    cancer_kmers, cancer_gene_hits, cancer_regions, gene_kmer_counts = find_cancer_kmers(genome, k, genes, cancer_genes)
    print(f"Found {len(cancer_kmers)} unique k-mers in cancer genes")
    top_cancer = [item for item in sorted(cancer_kmers.items(), key=lambda x: x[1], reverse=True) if item[1] < 10000][:5]
    print("Top 5 cancer gene k-mers:")
    for kmer, count in top_cancer:
        genes_hit = cancer_gene_hits[kmer]
        unique_genes = list(set(genes_hit))[:3]
        print(f"{kmer}: {count} (in genes: {', '.join(unique_genes)}{'...' if len(set(genes_hit)) > 3 else ''})")

    # Save results
    output_file = f"kmer_results_{args.filename.split('.')[0]}.txt"
    with open(output_file, "w") as f:
        f.write(f"File: {args.filename}\nBasic: {time_b:.4f} sec\nParallel: {time_p:.4f} sec\nSpeedup: {speedup:.2f}%\n")
        f.write(f"Cancer k-mers: {len(cancer_kmers)}\nTop 5:\n")
        for kmer, count in top_cancer:
            genes_hit = cancer_gene_hits[kmer]
            unique_genes = list(set(genes_hit))[:3]
            f.write(f"{kmer}: {count} (in genes: {', '.join(unique_genes)}{'...' if len(set(genes_hit)) > 3 else ''})\n")

    # Plot 1: Top Cancer Gene K-mers
    kmers, counts = zip(*top_cancer)
    plt.figure(figsize=(10, 6))
    plt.bar(kmers, counts, color='skyblue')
    plt.xlabel("K-mers")
    plt.ylabel("Frequency")
    plt.title("Top 5 Cancer Gene K-mers (Chr1)")
    plt.xticks(rotation=45)
    for i, count in enumerate(counts):
        plt.text(i, count + 200, f"{count}", ha='center', va='bottom')
    plt.tight_layout()
    plt.savefig(f"top_cancer_kmers_{args.filename.split('.')[0]}.png")
    plt.close()
    print(f"Saved top k-mers plot as 'top_cancer_kmers_{args.filename.split('.')[0]}.png'")

    # Plot 2: Runtime Comparison
    plt.figure(figsize=(8, 5))
    runtimes = [time_b, time_p]
    labels = ['Basic', 'Parallel']
    plt.bar(labels, runtimes, color=['gray', 'teal'])
    plt.ylabel("Runtime (seconds)")
    plt.title("Runtime Comparison (Chr1 Analysis)")
    for i, runtime in enumerate(runtimes):
        plt.text(i, runtime + 1, f"{runtime:.2f}s", ha='center', va='bottom')
    plt.text(0.5, max(runtimes) * 0.9, f"Speedup: {speedup:.2f}%", ha='center', va='center', bbox=dict(facecolor='white', alpha=0.5))
    plt.tight_layout()
    plt.savefig(f"runtime_comparison_{args.filename.split('.')[0]}.png")
    plt.close()
    print(f"Saved runtime plot as 'runtime_comparison_{args.filename.split('.')[0]}.png'")

    # Plot 3: Pie Chart of Cancer vs. Non-Cancer Bases
    total_bases = len(genome)
    cancer_bases = len(cancer_regions)
    non_cancer_bases = total_bases - cancer_bases
    plt.figure(figsize=(7, 7))
    plt.pie([cancer_bases, non_cancer_bases], labels=['Cancer Genes', 'Non-Cancer'], autopct='%1.1f%%', colors=['salmon', 'lightgray'])
    plt.title("Cancer vs. Non-Cancer Bases (Chr1)")
    plt.savefig(f"cancer_base_proportion_{args.filename.split('.')[0]}.png")
    plt.close()
    print(f"Saved pie chart as 'cancer_base_proportion_{args.filename.split('.')[0]}.png'")

    # Plot 4: Heatmap of K-mer Frequency vs. Gene
    # Select top 5 k-mers and top 5 genes by total k-mer count
    top_kmers = [kmer for kmer, _ in top_cancer]
    gene_totals = defaultdict(int)
    for gene in cancer_genes:
        for kmer in gene_kmer_counts[gene]:
            gene_totals[gene] += gene_kmer_counts[gene][kmer]
    top_genes = sorted(gene_totals.items(), key=lambda x: x[1], reverse=True)[:5]
    top_gene_names = [gene for gene, _ in top_genes]

    # Build heatmap data
    heatmap_data = np.zeros((len(top_kmers), len(top_gene_names)))
    for i, kmer in enumerate(top_kmers):
        for j, gene in enumerate(top_gene_names):
            heatmap_data[i, j] = gene_kmer_counts[gene].get(kmer, 0)

    # Create heatmap
    plt.figure(figsize=(12, 8))
    plt.imshow(heatmap_data, cmap='YlOrRd', aspect='auto')
    plt.colorbar(label='Frequency')
    plt.xticks(np.arange(len(top_gene_names)), top_gene_names, rotation=45, ha='right')
    plt.yticks(np.arange(len(top_kmers)), top_kmers)
    plt.xlabel("Cancer Genes")
    plt.ylabel("K-mers")
    plt.title("K-mer Frequency in Top Cancer Genes (Chr1)")
    plt.tight_layout()
    plt.savefig(f"kmer_gene_heatmap_{args.filename.split('.')[0]}.png")
    plt.close()
    print(f"Saved heatmap as 'kmer_gene_heatmap_{args.filename.split('.')[0]}.png'")

if __name__ == "__main__":
    main()
    
