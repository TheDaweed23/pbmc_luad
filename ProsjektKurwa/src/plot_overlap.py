# src/plot_overlap.py
from pathlib import Path
import math
import pandas as pd
import matplotlib.pyplot as plt

# Paths
ROOT = Path(__file__).resolve().parents[1]
OVER = ROOT / "results" / "overlap_pathways"
FIGS = OVER / "figures"
FIGS.mkdir(parents=True, exist_ok=True)

# Colors per library
PALETTE = {
    "GO_Biological_Process_2023": "#1f77b4",
    "KEGG_2021_Human": "#ff7f0e",
    "Reactome_2022": "#2ca02c",
}

def load_overlap(name: str) -> pd.DataFrame:
    """Load overlap CSV, add label, -log10(max p), and size."""
    path = OVER / f"{name}_pathway_overlap.csv"
    df = pd.read_csv(path)
    df["term"] = df["Term_norm"].str.replace(r"\s*\[.*\]$", "", regex=True)
    df["neglog10_maxP"] = df["max_adjP"].apply(lambda p: -math.log10(p) if p > 0 else 300)
    df["size"] = df["combined_score_sum"]
    return df

def top_per_library(df: pd.DataFrame, n=10) -> pd.DataFrame:
    """Top n per library by lowest max_adjP, then highest combined score."""
    df = df.sort_values(["max_adjP", "combined_score_sum"], ascending=[True, False])
    return df.groupby("library", as_index=False).head(n)

def dotplot(df: pd.DataFrame, title: str, outpng: Path):
    """Dot plot of top terms per library."""
    dft = top_per_library(df, n=10).copy()
    dft["rank"] = dft.groupby("library")["max_adjP"].rank(method="first")
    dft = dft.sort_values(["library", "rank"])

    y_labels = list(dict.fromkeys(dft["term"].tolist()))
    y_index = {t: i for i, t in enumerate(y_labels)}

    plt.figure(figsize=(10, 0.4 * len(y_labels) + 2))
    for _, r in dft.iterrows():
        y = y_index[r["term"]]
        color = PALETTE.get(r["library"], "#7f7f7f")
        s = max(20, min(300, r["size"]))  # cap size for readability
        plt.scatter(r["neglog10_maxP"], y, s=s, c=color, alpha=0.8, edgecolor="k", linewidths=0.3)

    plt.yticks(range(len(y_labels)), y_labels)
    plt.xlabel(r"$-\log_{10}(\mathrm{max\ adjusted\ }p)$")
    plt.title(title)

    # Legend (only for libraries present)
    handles, labels = [], []
    for lib, color in PALETTE.items():
        if lib in dft["library"].values:
            handles.append(plt.Line2D([], [], marker="o", linestyle="None", color=color, markeredgecolor="k", markersize=8))
            labels.append(lib.replace("_", " "))
    if handles:
        plt.legend(handles, labels, loc="lower right", frameon=True)

    plt.grid(axis="x", linestyle=":", alpha=0.4)
    plt.tight_layout()
    plt.savefig(outpng, dpi=300)
    plt.close()

def counts_by_library(dfa: pd.DataFrame, dfb: pd.DataFrame):
    """Bar plot of overlap counts per library for both directions."""
    a = dfa.groupby("library").size().rename("exUp_caDown")
    b = dfb.groupby("library").size().rename("exDown_caUp")
    g = pd.concat([a, b], axis=1).fillna(0).astype(int)

    ax = g.plot(kind="bar", figsize=(8, 4))
    ax.set_ylabel("overlapping pathways")
    ax.set_title("overlap counts by library")
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.savefig(FIGS / "overlap_counts_by_library.png", dpi=300)
    plt.close()

def main():
    exUp_caDown = load_overlap("exercise_up__cancer_down")
    exDown_caUp = load_overlap("exercise_down__cancer_up")

    dotplot(
        exUp_caDown,
        "overlapping pathways: exercise↑  ∩  cancer↓",
        FIGS / "exUp_caDown_dotplot.png",
    )
    dotplot(
        exDown_caUp,
        "overlapping pathways: exercise↓  ∩  cancer↑",
        FIGS / "exDown_caUp_dotplot.png",
    )
    counts_by_library(exUp_caDown, exDown_caUp)

if __name__ == "__main__":
    main()
