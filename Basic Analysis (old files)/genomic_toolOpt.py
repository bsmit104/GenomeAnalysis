# Brayden Smith
# jbrayden35@gmail.com
# Early script used to experiment with optimizing run times

import time
import random

def generate_mock_genome(length):
    bases = ['A', 'T', 'C', 'G']
    return ''.join(random.choice(bases) for _ in range(length))

genome = generate_mock_genome(1000000)
substring = "ATG"

def find_positions_basic(genome, substring):
    start_time = time.time()
    positions = []
    for i in range(len(genome) - len(substring) + 1):
        if genome[i:i + len(substring)] == substring:
            positions.append(i)
    end_time = time.time()
    return positions, end_time - start_time

# Test baseline
pos_basic, time_basic = find_positions_basic(genome, substring)
print(f"Basic: Found {len(pos_basic)} occurrences in {time_basic:.4f} seconds")


def find_positions_optimized(genome, substring):
    start_time = time.time()
    positions = []
    pos = -1
    while True:
        pos = genome.find(substring, pos + 1)
        if pos == -1:  # No more occurrences
            break
        positions.append(pos)
    end_time = time.time()
    return positions, end_time - start_time

# Test optimized
pos_opt, time_opt = find_positions_optimized(genome, substring)
print(f"Optimized: Found {len(pos_opt)} occurrences in {time_opt:.4f} seconds")


print("Running comparisons...")
pos_b, time_b = find_positions_basic(genome, substring)
pos_o, time_o = find_positions_optimized(genome, substring)

print(f"Basic: {time_b:.4f} sec, Count: {len(pos_b)}")
print(f"Optimized: {time_o:.4f} sec, Count: {len(pos_o)}")
speedup = (time_b - time_o) / time_b * 100
print(f"Speedup: {speedup:.2f}%")


# BASIC: Loops over every position and slices the string—slow 
# for large genomes (~0.1-0.2 sec expected for 1M bases).
# OPTIMIZED: Uses str.find(), which is implemented in C and 
# skips unnecessary checks by starting from the last found 
# position. This cuts runtime significantly (likely ~0.05-0.1 sec).
# SPEED UP: With 1M bases, you should see 40-60% improvement, 
# exceeding your 30% goal.