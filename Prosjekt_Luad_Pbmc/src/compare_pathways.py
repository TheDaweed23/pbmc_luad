# src/compare_pathways.py
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ENR = ROOT / "results" / "enrichment"
OUT = ROOT / "results" / "overlap_pathways"
OUT.mkdir(parents=True, exist_ok=True)


def load_sig(df_path: Path, q: float) -> pd.DataFrame:
    """Load enrichment CSV, keep rows with adj p ≤ q, add normalized 'Term_norm'."""
    df = pd.read_csv(df_path)
    df = df[df["Adjusted P-value"] <= q].copy()
    df["Term_norm"] = df["Term"].str.replace(r"\s*\(.*\)$", "", regex=True).str.strip()
    return df


def chosen_q(dir_: Path) -> float:
    """Return q from '*_chosen_q.txt' if present, else 0.05."""
    q_txt = sorted(dir_.glob("*_chosen_q.txt"))
    if not q_txt:
        return 0.05
    return float(q_txt[0].read_text().strip())


def load_set(set_dir: Path, label: str) -> pd.DataFrame:
    """Load '{label}_raw_enrichment.csv' using directory-specific q."""
    q = chosen_q(set_dir)
    raw = sorted(set_dir.glob(f"{label}_raw_enrichment.csv"))[0]
    df = load_sig(raw, q)
    df["source_set"] = label
    return df


def intersect(a: pd.DataFrame, b: pd.DataFrame, name: str) -> pd.DataFrame:
    """Overlap on ['library','Term_norm'], compute simple combined stats, save CSV."""
    key = ["library", "Term_norm"]
    common = a.merge(b, on=key, suffixes=("_A", "_B"))
    common["max_adjP"] = common[["Adjusted P-value_A", "Adjusted P-value_B"]].max(axis=1)
    common["combined_score_sum"] = common["Combined Score_A"] + common["Combined Score_B"]

    cols = [
        "library", "Term_norm",
        "Adjusted P-value_A", "Adjusted P-value_B",
        "Odds Ratio_A", "Odds Ratio_B",
        "Combined Score_A", "Combined Score_B",
        "max_adjP", "combined_score_sum",
        "Genes_A", "Genes_B",
    ]
    common = common[cols].sort_values(["library", "max_adjP", "combined_score_sum"])
    common.to_csv(OUT / f"{name}_pathway_overlap.csv", index=False)
    return common


def main():
    pbmc_dir = ENR / "pbmc"
    luad_dir = ENR / "luad"

    pbmc_up = load_set(pbmc_dir, "PBMC_up")
    pbmc_down = load_set(pbmc_dir, "PBMC_down")
    luad_up = load_set(luad_dir, "LUAD_up_targets")
    luad_down = load_set(luad_dir, "LUAD_down_targets")

    ex_up_ca_down = intersect(pbmc_up, luad_down, "exercise_up__cancer_down")
    ex_down_ca_up = intersect(pbmc_down, luad_up, "exercise_down__cancer_up")

    # top 20 per library quick lists
    for df, tag in [(ex_up_ca_down, "exUp_caDown"), (ex_down_ca_up, "exDown_caUp")]:
        (OUT / f"TOP20_{tag}.csv").write_text(
            df.groupby("library").head(20).to_csv(index=False)
        )


if __name__ == "__main__":
    main()
