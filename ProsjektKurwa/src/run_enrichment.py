# src/run_enrichment.py
from pathlib import Path
import pandas as pd
from gseapy import enrichr

# Paths/config
PROJECT_ROOT = Path(__file__).resolve().parents[1]
GENE_DIR = PROJECT_ROOT / "filtered_gene_lists"
COMP_DIR = PROJECT_ROOT / "results" / "enrichment" / "comparison"
OUT_DIR = PROJECT_ROOT / "results" / "enrichment"
(OUT_DIR / "pbmc").mkdir(parents=True, exist_ok=True)
(OUT_DIR / "luad").mkdir(parents=True, exist_ok=True)
COMP_DIR.mkdir(parents=True, exist_ok=True)

LIBRARIES = [
    "KEGG_2021_Human",
    "Reactome_2022",
    "GO_Biological_Process_2023",
]
ORGANISM = "Human"


def read_gene_list(path: Path) -> list[str]:
    """Read non-empty, stripped lines; de-duplicate preserving order."""
    genes = []
    with path.open() as f:
        for line in f:
            g = line.strip()
            if g:
                genes.append(g)
    return list(dict.fromkeys(genes))


def run_one_enrich(gene_list: list[str], label: str, outdir: Path) -> pd.DataFrame:
    """Run Enrichr across LIBRARIES, tag results, save raw CSV."""
    outdir.mkdir(parents=True, exist_ok=True)
    all_res = []
    for lib in LIBRARIES:
        res = enrichr(
            gene_list=gene_list,
            gene_sets=lib,
            organism=ORGANISM,
            outdir=None,   # manual saving
            cutoff=1.0,    # no filtering; filter later
        )
        if res is None or res.results is None or res.results.empty:
            continue
        df = res.results.copy()
        df["library"] = lib
        df["set_label"] = label
        all_res.append(df)

    if not all_res:
        return pd.DataFrame()

    df = pd.concat(all_res, ignore_index=True)
    keep = [
        "Term", "Adjusted P-value", "P-value", "Odds Ratio", "Combined Score",
        "Overlap", "Genes", "library", "set_label",
    ]
    df = df[keep]
    df.to_csv(outdir / f"{label}_raw_enrichment.csv", index=False)
    return df


def auto_threshold(df: pd.DataFrame, target_min=10, target_max=200) -> float:
    """
    Choose q (Adjusted P-value cutoff) to yield a result count in [min, max].
    """
    if df.empty:
        return 0.05
    candidates = [0.01, 0.02, 0.03, 0.05, 0.07, 0.1]
    for q in candidates:
        n = (df["Adjusted P-value"] <= q).sum()
        if target_min <= n <= target_max:
            return q
    if (df["Adjusted P-value"] <= 0.05).sum() > 0:
        return 0.05
    return 0.1


def save_filtered(df: pd.DataFrame, q: float, outdir: Path, label: str):
    """Filter by q, sort, save filtered CSV and chosen_q marker."""
    f = df[df["Adjusted P-value"] <= q].copy()
    f.sort_values(["Adjusted P-value", "library"], inplace=True)
    f.to_csv(outdir / f"{label}_enrichment_q{q:.2f}.csv", index=False)
    (outdir / f"{label}_chosen_q.txt").write_text(str(q))


def main():
    inputs = {
        "PBMC_up": GENE_DIR / "pbmc_up.txt",
        "PBMC_down": GENE_DIR / "pbmc_down.txt",
        "LUAD_up_targets": COMP_DIR / "luad_up_targets.txt",
        "LUAD_down_targets": COMP_DIR / "luad_down_targets.txt",
    }

    for label, path in inputs.items():
        genes = read_gene_list(path)
        if not genes:
            print(f"[warn] {label}: input list is empty at {path}")
            continue

        out_bucket = OUT_DIR / ("pbmc" if label.startswith("PBMC") else "luad")
        raw = run_one_enrich(genes, label, out_bucket)
        if raw.empty:
            print(f"[warn] {label}: enrichment returned no terms")
            continue

        q = auto_threshold(raw)
        save_filtered(raw, q, out_bucket, label)
        print(f"[ok] {label}: saved filtered enrichment with q={q:.2f} to {out_bucket}")


if __name__ == "__main__":
    main()
