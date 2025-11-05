# src/overlap_genes.py
from pathlib import Path

def read_list(p: Path) -> set[str]:
    """Read a text file into a set of non-empty, stripped lines."""
    return set(x.strip() for x in p.read_text().splitlines() if x.strip())


# Define paths
ROOT = Path(__file__).resolve().parents[1]
GL = ROOT / "filtered_gene_lists"
COMP = ROOT / "results" / "enrichment" / "comparison"
OUT = ROOT / "results" / "summary_tables"
OUT.mkdir(parents=True, exist_ok=True)

# Load gene lists
pbmc_up = read_list(GL / "pbmc_up.txt")
pbmc_down = read_list(GL / "pbmc_down.txt")
luad_up_t = read_list(COMP / "luad_up_targets.txt")
luad_down_t = read_list(COMP / "luad_down_targets.txt")

# Find overlaps
exUp_caDown = sorted(pbmc_up & luad_down_t)
exDown_caUp = sorted(pbmc_down & luad_up_t)

# Save results
(OUT / "gene_overlap_exUp_caDown.txt").write_text("\n".join(exUp_caDown))
(OUT / "gene_overlap_exDown_caUp.txt").write_text("\n".join(exDown_caUp))

# Print summary
print(f"exercise↑ ∩ cancer↓: {len(exUp_caDown)} genes -> results/summary_tables/gene_overlap_exUp_caDown.txt")
print(f"exercise↓ ∩ cancer↑: {len(exDown_caUp)} genes -> results/summary_tables/gene_overlap_exDown_caUp.txt")
