# src/plot_gene_overlap_counts.py
from pathlib import Path
import matplotlib.pyplot as plt

# Paths
ROOT = Path(__file__).resolve().parents[1]
SUM = ROOT / "results" / "summary_tables"
FIGS = SUM / "figures"
FIGS.mkdir(parents=True, exist_ok=True)

a_path = SUM / "gene_overlap_exUp_caDown.txt"
b_path = SUM / "gene_overlap_exDown_caUp.txt"

def read_list(p: Path) -> list[str]:
    """Return non-empty, stripped lines from a text file."""
    return [x.strip() for x in p.read_text(encoding="utf-8").splitlines() if x.strip()]

def main():
    # Load overlaps
    a = read_list(a_path)
    b = read_list(b_path)

    # Bar plot: counts of overlaps
    bar_labels = ["exercise↑ ∩ cancer↓", "exercise↓ ∩ cancer↑"]
    counts = [len(a), len(b)]

    plt.figure(figsize=(6, 4))
    plt.bar(bar_labels, counts)
    for i, c in enumerate(counts):
        plt.text(i, c + max(counts) * 0.02, str(c), ha="center", va="bottom")
    plt.ylabel("overlap gene count")
    plt.title("gene-level overlaps")
    plt.tight_layout()
    plt.savefig(FIGS / "gene_overlap_counts.png", dpi=300)
    plt.close()

    # ASCII-only text preview (first 30 genes per group)
    preview = []
    preview.append("exercise_up AND cancer_down (first 30):\n" + ", ".join(a[:30]) + "\n")
    preview.append("exercise_down AND cancer_up (first 30):\n" + ", ".join(b[:30]) + "\n")
    (FIGS / "gene_overlap_lists_preview.txt").write_text("\n".join(preview), encoding="utf-8")

if __name__ == "__main__":
    main()
