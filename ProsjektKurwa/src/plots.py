from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re

# -------------------------------
ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "results" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Point directly to the Enrichr CSVs
EXERCISE_UP   = ROOT / "results" / "enrichment" / "pbmc" / "PBMC_up_enrichment_q0.10.csv"
EXERCISE_DOWN = ROOT / "results" / "enrichment" / "pbmc" / "PBMC_down_enrichment_q0.10.csv"
CANCER_UP     = ROOT / "results" / "enrichment" / "luad" / "LUAD_up_targets_enrichment_q0.05.csv"
CANCER_DOWN   = ROOT / "results" / "enrichment" / "luad" / "LUAD_down_targets_enrichment_q0.05.csv"

# Optional fallbacks if the q-threshold files don’t exist
EXERCISE_UP_FALLBACK   = ROOT / "results" / "enrichment" / "pbmc" / "PBMC_up_raw_enrichment.csv"
EXERCISE_DOWN_FALLBACK = ROOT / "results" / "enrichment" / "pbmc" / "PBMC_down_raw_enrichment.csv"
CANCER_UP_FALLBACK     = ROOT / "results" / "enrichment" / "luad" / "LUAD_up_targets_raw_enrichment.csv"
CANCER_DOWN_FALLBACK   = ROOT / "results" / "enrichment" / "luad" / "LUAD_down_targets_raw_enrichment.csv"

# Themes (regex, case-insensitive)
THEMES = {
    "immune": [
        r"\bimmune\b", r"\bimmunity\b",
        r"T[- ]?cell|B[- ]?cell|lymphocyte|TCR|CD[0-9]+",
        r"NK|natural killer|cytotox|granzyme|perforin",
        r"antigen|MHC|HLA|presentation",
        r"interferon|IFN[- ]?(alpha|beta|gamma)|STAT1|IRF"
    ],
    "inflammatory": [
        r"IL-6|JAK|STAT3|TNF|NF[- ]?KB|S100A8|S100A9",
        r"acute inflammatory|inflammation|chemokine|complement|coagulation"
    ],
    "metabolic": [
        r"oxidative phosphorylation|OXPHOS|electron transport",
        r"fatty[- ]?acid|beta[- ]?oxidation|PPAR|mitochond",
        r"glycolysis|TCA|tricarboxylic"
    ],
}

# Immune-only filter for the dotplot
IMMUNE_REGEX = re.compile("|".join(THEMES["immune"]), re.IGNORECASE)


def _first_existing(*paths: Path) -> Path | None:
    for p in paths:
        if p and p.exists():
            return p
    return None


def _load_enrichr(path: Path) -> pd.DataFrame:
    """Load Enrichr-like CSV robustly; returns columns: Term, adjp, Overlap, mlog10_adj."""
    if path is None or not path.exists():
        return pd.DataFrame(columns=["Term", "adjp", "Overlap", "mlog10_adj"])

    df = pd.read_csv(path)
    cols = {c.lower().strip(): c for c in df.columns}

    # Term column
    term_col = None
    for k in ["term", "terms", "name", "pathway", "description"]:
        if k in cols:
            term_col = cols[k]; break
    if term_col is None:
        raise ValueError(f"Could not find a 'Term' column in {path.name}")

    # Adjusted p-value (or fallback to raw p)
    adjp_col = None
    for k in ["adjusted p-value", "adj p", "adjp", "fdr", "q-value", "qvalue", "p_adj"]:
        if k in cols:
            adjp_col = cols[k]; break
    if adjp_col is None:
        for k in ["p-value", "pvalue", "p value"]:
            if k in cols:
                adjp_col = cols[k]; break
    if adjp_col is None:
        raise ValueError(f"Could not find a p-value column in {path.name}")

    # Overlap / gene count (optional)
    overlap_col = None
    for k in ["overlap", "overlapping genes", "genes", "gene count", "count"]:
        if k in cols:
            overlap_col = cols[k]; break

    out = pd.DataFrame({
        "Term": df[term_col].astype(str),
        "adjp": pd.to_numeric(df[adjp_col], errors="coerce")
    })
    if overlap_col:
        ov = df[overlap_col].astype(str)
        num = ov.str.extract(r"(\d+)", expand=False).astype(float)
        out["Overlap"] = num
    else:
        out["Overlap"] = np.nan

    out = out.dropna(subset=["adjp"]).replace(0, np.nextafter(0, 1))
    out["mlog10_adj"] = -np.log10(out["adjp"])
    return out


def _theme_score(df: pd.DataFrame, regex_list: list[str]) -> float:
    """Median −log10(adjp) for terms matching any regex (0.0 if none)."""
    if df is None or df.empty:
        return 0.0
    mask = np.zeros(len(df), dtype=bool)
    for rx in regex_list:
        mask |= df["Term"].str.contains(rx, regex=True, case=False, na=False)
    return float(df.loc[mask, "mlog10_adj"].median()) if mask.any() else 0.0


def _clean_term(t: str) -> str:
    # strip Reactome IDs like "R-HSA-123456" and extra spaces
    t = re.sub(r"\sR-HSA-\d+\b", "", t)
    t = re.sub(r"\s\(GO:\d+\)", "", t)
    return t.strip()

def plot_dotplot(ex_up: pd.DataFrame, ca_dn: pd.DataFrame, save_path: Path, top_n=12):
    # broaden “immune-like” selection; if empty, fall back to top_n by adjp
    ex_df = ex_up[ex_up["Term"].str.contains(IMMUNE_REGEX)]
    ca_df = ca_dn[ca_dn["Term"].str.contains(IMMUNE_REGEX)]
    if ex_df.empty: ex_df = ex_up.sort_values("adjp").head(top_n).copy()
    if ca_df.empty: ca_df = ca_dn.sort_values("adjp").head(top_n).copy()

    ex_df = ex_df.sort_values("adjp").head(top_n).copy()
    ca_df = ca_df.sort_values("adjp").head(top_n).copy()
    ex_df["Set"] = "Exercise↑"; ca_df["Set"] = "Cancer↓"

    comb = pd.concat([ex_df, ca_df], ignore_index=True)
    comb["Term"] = comb["Term"].map(_clean_term)

    # order y by average rank across sets so similar terms sit near each other
    comb["rank"] = comb.groupby("Set")["adjp"].rank(method="first")
    order = (comb.groupby("Term")["rank"].mean().sort_values()).index.tolist()
    y_map = {t: i for i, t in enumerate(order)}

    plt.figure(figsize=(10, 7))
    for s, color, marker in [("Exercise↑", "tab:blue", "o"), ("Cancer↓", "tab:orange", "s")]:
        sub = comb[comb["Set"] == s]
        x = sub["mlog10_adj"].values
        y = [y_map[t] + (0.15 if s == "Cancer↓" else -0.15) for t in sub["Term"]]
        size = (sub["Overlap"].fillna(sub["Overlap"].median() if not sub["Overlap"].dropna().empty else 5) + 2) * 8
        plt.scatter(x, y, s=size, c=color, marker=marker, alpha=0.9, edgecolors="k", linewidths=0.3, label=s)

    plt.yticks(range(len(order)), order)
    plt.xlabel("−log10(adjusted p)")
    plt.title("Immune-related pathways: Exercise↑ vs Cancer↓ (ORA)")
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()
    print(f"Saved: {save_path}")



def plot_heatmatrix(ex_up, ex_dn, ca_up, ca_dn, save_path_png: Path, save_tbl_csv: Path):
    # theme × (Exercise↑, Exercise↓, Cancer↑, Cancer↓) with median −log10(adjp)
    rows = []
    for theme, regs in THEMES.items():
        rows.append({
            "Theme": theme,
            "Exercise↑": _theme_score(ex_up, regs),
            "Exercise↓": _theme_score(ex_dn, regs),
            "Cancer↑":   _theme_score(ca_up, regs),
            "Cancer↓":   _theme_score(ca_dn, regs),
        })
    mat = pd.DataFrame(rows).set_index("Theme")

    # Heatmap with annotations
    plt.figure(figsize=(8.2, 4.8))
    data = mat[["Exercise↑", "Exercise↓", "Cancer↑", "Cancer↓"]].values
    im = plt.imshow(data, aspect="auto", cmap="Blues")
    plt.colorbar(im, fraction=0.046, pad=0.04, label="Median −log10(adjp)")
    plt.xticks(range(4), ["Exercise↑", "Exercise↓", "Cancer↑", "Cancer↓"])
    plt.yticks(range(mat.shape[0]), mat.index)
    for (i, j), val in np.ndenumerate(data):
        plt.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=8, color="black")
    plt.title("Pathway-theme enrichment across Exercise/LUAD (ORA)")
    plt.tight_layout()
    plt.savefig(save_path_png, dpi=300)
    plt.close()
    print(f"Saved: {save_path_png}")

    # Direction summary useful for table1
    dir_score = (mat["Exercise↑"] + mat["Cancer↓"]) - (mat["Exercise↓"] + mat["Cancer↑"])
    if dir_score.abs().max() > 0:
        dir_norm = dir_score / dir_score.abs().max()
    else:
        dir_norm = dir_score
    pd.DataFrame({"Theme": mat.index, "DirectionScore": dir_norm.values}) \
      .to_csv(save_tbl_csv, index=False)
    print(f"Saved: {save_tbl_csv}")

def _clean_term(t: str) -> str:
    # strip Reactome/GO IDs and extra spaces
    t = re.sub(r"\sR-HSA-\d+\b", "", t)
    t = re.sub(r"\s\(GO:\d+\)", "", t)
    t = re.sub(r"\s+\bHomo sapiens\b", "", t, flags=re.I)
    return t.strip()

def plot_mirror_bars(ex_up: pd.DataFrame, ca_dn: pd.DataFrame, save_path: Path, top_n=12):
    """
    Mirrored horizontal bar chart:
      - left (negative x): Exercise↑  (blue)
      - right (positive x): Cancer↓   (orange)
    x-axis is −log10(adjusted p). Same terms aligned on y.
    """
    # prefer immune-like terms; fall back to top by significance if empty
    ex_df = ex_up[ex_up["Term"].str.contains(IMMUNE_REGEX, regex=True, na=False)]
    ca_df = ca_dn[ca_dn["Term"].str.contains(IMMUNE_REGEX, regex=True, na=False)]
    if ex_df.empty: ex_df = ex_up.sort_values("adjp").head(top_n).copy()
    if ca_df.empty: ca_df = ca_dn.sort_values("adjp").head(top_n).copy()

    # keep top_n most significant per side
    ex_df = ex_df.sort_values("adjp").head(top_n).copy()
    ca_df = ca_df.sort_values("adjp").head(top_n).copy()

    ex_df["Term"] = ex_df["Term"].map(_clean_term)
    ca_df["Term"] = ca_df["Term"].map(_clean_term)

    # union of terms, ordered by combined rank so similar items sit together
    ex_df["rank_ex"] = ex_df["adjp"].rank(method="first")
    ca_df["rank_ca"] = ca_df["adjp"].rank(method="first")
    joined = pd.merge(ex_df[["Term","mlog10_adj","rank_ex"]].rename(columns={"mlog10_adj":"ex_mlog"}),
                      ca_df[["Term","mlog10_adj","rank_ca"]].rename(columns={"mlog10_adj":"ca_mlog"}),
                      on="Term", how="outer")
    joined["rank_comb"] = joined[["rank_ex","rank_ca"]].min(axis=1)
    joined = joined.sort_values(["rank_comb","Term"]).head(top_n)

    terms = joined["Term"].tolist()
    ex_vals = -joined["ex_mlog"].fillna(0).values   # negative for left bars
    ca_vals =  joined["ca_mlog"].fillna(0).values   # positive for right bars

    # figure
    import matplotlib.pyplot as plt
    plt.figure(figsize=(11, 7))
    y = np.arange(len(terms))

    plt.barh(y, ex_vals, color="tab:blue", alpha=0.9, label="Exercise↑")
    plt.barh(y, ca_vals, color="tab:orange", alpha=0.9, label="Cancer↓")

    # labels at bar ends
    def _annot(vals, offset=0.15):
        for yi, v in enumerate(vals):
            if v == 0: continue
            plt.text(v + (offset if v>0 else -offset), yi, f"{abs(v):.2f}",
                     va="center", ha="left" if v>0 else "right", fontsize=8)

    _annot(ex_vals); _annot(ca_vals)

    plt.yticks(y, terms)
    plt.axvline(0, color="k", linewidth=0.8)
    # symmetric x-limits
    xmax = max(abs(ex_vals).max(), abs(ca_vals).max()) * 1.15
    plt.xlim(-xmax, xmax)
    plt.xlabel("−log10(adjusted p)")
    plt.title("Immune-related pathways (ORA): Exercise↑ vs Cancer↓")
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    plt.close()
    print(f"Saved: {save_path}")

def main():
    # Resolve with fallbacks
    paths = {
        "exercise_up":   _first_existing(EXERCISE_UP, EXERCISE_UP_FALLBACK),
        "exercise_down": _first_existing(EXERCISE_DOWN, EXERCISE_DOWN_FALLBACK),
        "cancer_up":     _first_existing(CANCER_UP, CANCER_UP_FALLBACK),
        "cancer_down":   _first_existing(CANCER_DOWN, CANCER_DOWN_FALLBACK),
    }
    print("Using enrichment files:")
    for k, v in paths.items():
        print(f"  {k}: {v}")

    ex_up = _load_enrichr(paths["exercise_up"])
    ex_dn = _load_enrichr(paths["exercise_down"])
    ca_up = _load_enrichr(paths["cancer_up"])
    ca_dn = _load_enrichr(paths["cancer_down"])

    # 1) Immune figure: mirrored bars
    if not ex_up.empty and not ca_dn.empty:
        plot_mirror_bars(ex_up, ca_dn, FIG_DIR / "pathway_mirror_bars.png", top_n=12)
    else:
        print("Mirror bars skipped: missing Exercise↑ or Cancer↓ terms.")

    # 2) Theme heatmatrix
    plot_heatmatrix(ex_up, ex_dn, ca_up, ca_dn,
                    FIG_DIR / "pathway_heatmatrix.png",
                    FIG_DIR / "pathway_direction_summary.csv")


if __name__ == "__main__":
    main()
