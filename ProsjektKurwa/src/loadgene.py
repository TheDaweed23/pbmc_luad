import os

# Base paths (relative to /src/)
base_dir = os.path.dirname(__file__)
gene_list_dir = os.path.join(base_dir, "..", "filtered_gene_lists")

# Input files
pbmc_up_file = os.path.join(gene_list_dir, "pbmc_up.txt")
pbmc_down_file = os.path.join(gene_list_dir, "pbmc_down.txt")
luad_up_file = os.path.join(gene_list_dir, "luad_up.txt")
luad_down_file = os.path.join(gene_list_dir, "luad_down.txt")

# Load PBMC gene lists
with open(pbmc_up_file, "r") as f:
    pbmc_up_genes = [line.strip() for line in f if line.strip()]

with open(pbmc_down_file, "r") as f:
    pbmc_down_genes = [line.strip() for line in f if line.strip()]

# Load LUAD miRNA lists
with open(luad_up_file, "r") as f:
    luad_up_mirs = [line.strip() for line in f if line.strip()]

with open(luad_down_file, "r") as f:
    luad_down_mirs = [line.strip() for line in f if line.strip()]

# Print summary
print(f"Loaded {len(pbmc_up_genes)} PBMC-up genes")
print(f"Loaded {len(pbmc_down_genes)} PBMC-down genes")
print(f"Loaded {len(luad_up_mirs)} LUAD-up miRNAs")
print(f"Loaded {len(luad_down_mirs)} LUAD-down miRNAs")
