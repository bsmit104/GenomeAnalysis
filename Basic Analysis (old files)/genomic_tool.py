import random

def generate_mock_genome(length):
    bases = ['A', 'T', 'C', 'G']
    return ''.join(random.choice(bases) for _ in range(length))

# Fake sample sequence
genome = generate_mock_genome(10000)
print(genome[:50])



# Using Python’s built-in count() method as the “unoptimized” baseline.
import time

def count_substring_basic(genome, substring):
    start_time = time.time()
    count = genome.count(substring)
    end_time = time.time()
    return count, end_time - start_time

# Test 1
substring = "ATG"
count, runtime = count_substring_basic(genome, substring)
print(f"Basic method: Found {count} occurrences in {runtime:.4f} seconds")


genome = generate_mock_genome(1000000)
count, runtime = count_substring_basic(genome, substring)
print(f"Basic method: Found {count} occurrences in {runtime:.4f} seconds")
# Result : 3-5 times and take the average runtime (e.g., 0.015 seconds). This is the baseline.


def count_substring_optimized(genome, substring):
    start_time = time.time()
    count = 0
    sub_len = len(substring)
    for i in range(len(genome) - sub_len + 1):
        if genome[i:i + sub_len] == substring:
            count += 1
    end_time = time.time()
    return count, end_time - start_time

# Test 2
count_opt, runtime_opt = count_substring_optimized(genome, substring)
print(f"Optimized method: Found {count_opt} occurrences in {runtime_opt:.4f} seconds")


import numpy as np

def count_substring_numpy(genome, substring):
    start_time = time.time()
    # Convert genome to numpy array of characters
    genome_array = np.array(list(genome))
    sub_len = len(substring)
    # Use vectorized operations to check windows
    windows = [genome[i:i + sub_len] for i in range(len(genome) - sub_len + 1)]
    count = np.sum(np.array(windows) == substring)
    end_time = time.time()
    return count, end_time - start_time

# Test 3
count_np, runtime_np = count_substring_numpy(genome, substring)
print(f"NumPy method: Found {count_np} occurrences in {runtime_np:.4f} seconds")



print("Running comparisons...")
count_b, time_b = count_substring_basic(genome, substring)
count_o, time_o = count_substring_optimized(genome, substring)
count_n, time_n = count_substring_numpy(genome, substring)

print(f"Basic: {time_b:.4f} sec, Count: {count_b}")
print(f"Optimized: {time_o:.4f} sec, Count: {count_o}")
print(f"NumPy: {time_n:.4f} sec, Count: {count_n}")

speedup = (time_b - time_n) / time_b * 100
print(f"Speedup from Basic to NumPy: {speedup:.2f}%")



# OUTPUTS
with open("results.txt", "w") as f:
    f.write(f"Basic: {time_b:.4f} sec\nNumPy: {time_n:.4f} sec\nSpeedup: {speedup:.2f}%")


import matplotlib.pyplot as plt

methods = ['Basic', 'NumPy']
times = [time_b, time_n]
plt.bar(methods, times)
plt.ylabel("Runtime (seconds)")
plt.title("Genomic Analysis Runtime Comparison")
plt.savefig("runtime_comparison.png")

# ANALYSIS
# BASIC: 0.0030 sec: Fast because count() is efficient.
# OPTIMIZED: 0.1020 sec: 34x slower due to Python loop and slicing overhead.
# NUMPY: 0.4290 sec: 143x slower because it’s doing the 
# heavy lifting in Python first, then using NumPy inefficiently.
# SPEEDUP: -13984.50%: Negative because NumPy is slower, not faster. 
# (The formula (basic - numpy) / basic * 100 gives a huge negative value
# when numpy > basic.)

