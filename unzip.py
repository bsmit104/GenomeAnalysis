# Script for unzipping gz files for convenience

import gzip
import shutil

input_file = "Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz"
output_file = "human_chrom1.fa"

with gzip.open(input_file, "rb") as f_in:
    with open(output_file, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

print("File extracted successfully!")