import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier

from ai_config import RESULTS, RANDOM_STATE, TOP_K
from prep_luad_features import load_luad_data

def heatmap_topk():
    tables = RESULTS / "tables"
    plots = RESULTS / "plots"
    plots.mkdir(parents=True, exist_ok=True)

    X, y, _ = load_luad_data()
    # pick the same top-k used in training (if file exists)
    topk_file = tables / "topk_feature_scores.csv"
    if topk_file.exists():
        feats = pd.read_csv(topk_file)["feature"].head(TOP_K).tolist()
    else:
        feats = X.var().sort_values(ascending=False).head(TOP_K).index.tolist()

    Xs = X[feats]
    # z-score for heatmap
    Xz = pd.DataFrame(StandardScaler().fit_transform(Xs), index=Xs.index, columns=Xs.columns)
    # order samples by label
    order = y.sort_values().index
    Xz = Xz.loc[order]

    plt.figure(figsize=(8, 10))
    plt.imshow(Xz.values, aspect="auto", interpolation="nearest", cmap="coolwarm")
    plt.colorbar()
    plt.yticks(range(len(order)), order, fontsize=6)
    plt.xticks(range(len(feats)), feats, rotation=90, fontsize=6)
    plt.title("Top-K miRNA z-score heatmap (rows=samples, cols=features)")
    plt.tight_layout()
    plt.savefig(plots / "heatmap_topk.png", dpi=250)
    plt.close()

def pca_scatter():
    plots = RESULTS / "plots"
    X, y, _ = load_luad_data()
    Xz = StandardScaler().fit_transform(X)
    pca = PCA(n_components=2, random_state=RANDOM_STATE)
    Z = pca.fit_transform(Xz)
    plt.figure()
    for label, marker in [(0, "o"), (1, "s")]:
        pts = Z[y.values == label]
        plt.scatter(pts[:,0], pts[:,1], label=f"{'HC' if label==0 else 'LUAD'}", marker=marker)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend()
    plt.title("PCA of samples (all miRNAs)")
    plt.tight_layout()
    plt.savefig(plots / "pca_scatter.png", dpi=200)
    plt.close()

def rf_feature_importance():
    tables = RESULTS / "tables"
    plots = RESULTS / "plots"
    X, y, _ = load_luad_data()
    # use top-k features if available
    if (tables / "topk_feature_scores.csv").exists():
        feats = pd.read_csv(tables / "topk_feature_scores.csv")["feature"].head(TOP_K).tolist()
    else:
        feats = X.var().sort_values(ascending=False).head(TOP_K).index.tolist()
    Xs = X[feats]
    rf = RandomForestClassifier(n_estimators=600, random_state=0)
    rf.fit(Xs, y)
    imp = pd.Series(rf.feature_importances_, index=feats).sort_values(ascending=False)
    imp.head(30).to_csv(tables / "rf_feature_importance_top30.csv")
    plt.figure(figsize=(8, 6))
    imp.head(20).iloc[::-1].plot(kind="barh")
    plt.title("Random Forest feature importance (top 20)")
    plt.tight_layout()
    plt.savefig(plots / "rf_feature_importance_top20.png", dpi=220)
    plt.close()

if __name__ == "__main__":
    heatmap_topk()
    pca_scatter()
    rf_feature_importance()
