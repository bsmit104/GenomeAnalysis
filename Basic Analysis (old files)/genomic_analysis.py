import time
import random

# Generate mock genomic data
def generate_mock_genome(length):
    bases = ['A', 'T', 'C', 'G']
    return ''.join(random.choice(bases) for _ in range(length))

# Basic method: Naive loop with slicing
def find_positions_basic(genome, substring):
    start_time = time.time()
    positions = []
    for i in range(len(genome) - len(substring) + 1):
        if genome[i:i + len(substring)] == substring:
            positions.append(i)
    end_time = time.time()
    return positions, end_time - start_time

# Optimized method: Using str.find()
def find_positions_optimized(genome, substring):
    start_time = time.time()
    positions = []
    pos = -1
    while True:
        pos = genome.find(substring, pos + 1)
        if pos == -1:
            break
        positions.append(pos)
    end_time = time.time()
    return positions, end_time - start_time

# Main execution
def main():
    # Generate 1M base genome
    genome = generate_mock_genome(1000000)
    substring = "ATG"
    
    # Run comparisons
    print("Running comparisons...")
    pos_basic, time_basic = find_positions_basic(genome, substring)
    pos_opt, time_opt = find_positions_optimized(genome, substring)
    
    # Print results
    print(f"Basic: {time_basic:.4f} sec, Count: {len(pos_basic)}")
    print(f"Optimized: {time_opt:.4f} sec, Count: {len(pos_opt)}")
    speedup = (time_basic - time_opt) / time_basic * 100
    print(f"Speedup: {speedup:.2f}%")
    
    # Save results to file
    with open("genomic_results.txt", "w") as f:
        f.write(f"Basic: {time_basic:.4f} sec, Count: {len(pos_basic)}\n")
        f.write(f"Optimized: {time_opt:.4f} sec, Count: {len(pos_opt)}\n")
        f.write(f"Speedup: {speedup:.2f}%\n")

if __name__ == "__main__":
    main()