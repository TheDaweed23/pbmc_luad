import numpy as np
import pandas as pd
from typing import Dict, Tuple
from dataclasses import dataclass
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (accuracy_score, f1_score, roc_auc_score,
                             confusion_matrix, RocCurveDisplay)
import matplotlib.pyplot as plt
from pathlib import Path

from ai_config import RESULTS, TOP_K, TEST_SIZE, RANDOM_STATE, N_FOLDS, MODELS
from prep_luad_features import load_luad_data, read_miRNA_lists

@dataclass
class RunBundle:
    name: str                 # "topk" or "bio"
    X_train: pd.DataFrame
    X_test: pd.DataFrame
    y_train: pd.Series
    y_test: pd.Series
    selected_features: list

def _safe_folds(y, desired=5):
    # number of samples per class
    cls_counts = pd.Series(y).value_counts()
    min_cls = int(cls_counts.min())
    return max(2, min(desired, min_cls))

def _build_pipelines(n_features: int) -> Dict[str, Pipeline]:
    # scaler only where it matters
    logreg = Pipeline([
        ("scaler", StandardScaler(with_mean=True, with_std=True)),
        ("clf", LogisticRegression(max_iter=5000, solver="liblinear"))
    ])
    knn = Pipeline([
        ("scaler", StandardScaler(with_mean=True, with_std=True)),
        ("clf", KNeighborsClassifier(n_neighbors=5))
    ])
    rf = Pipeline([
        ("clf", RandomForestClassifier(n_estimators=400, random_state=RANDOM_STATE))
    ])
    out = {}
    if "logreg" in MODELS: out["logreg"] = logreg
    if "knn"   in MODELS: out["knn"]   = knn
    if "rf"    in MODELS: out["rf"]    = rf
    return out

def _split(X, y) -> Tuple[pd.DataFrame, pd.DataFrame, pd.Series, pd.Series]:
    return train_test_split(X, y, test_size=TEST_SIZE, random_state=RANDOM_STATE, stratify=y)

def _select_topk(X, y, k=TOP_K):
    # ensure finite
    Xc = X.replace([np.inf, -np.inf], np.nan).fillna(0.0)
    # drop any freshly zero-variance features
    var = Xc.var(axis=0)
    Xc = Xc.loc[:, var > 0]
    if Xc.shape[1] == 0:
        raise ValueError("No usable features after cleaning before ANOVA.")
    k = min(k, Xc.shape[1])
    selector = SelectKBest(score_func=f_classif, k=k)
    selector.fit(Xc, y)
    feats = Xc.columns[selector.get_support()].tolist()
    scores = selector.scores_
    pvals = selector.pvalues_
    return Xc[feats], feats, scores, pvals

def _select_bio(X, candidates):
    feats = [f for f in X.columns if f in set(candidates)]
    # keep at least something
    if len(feats) == 0:
        # fall back to all
        feats = X.columns.tolist()
    return X[feats], feats

def evaluate_and_plot(name: str, model_name: str, pipe: Pipeline,
                      X_train, y_train, X_test, y_test, outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    # Cross-val on train
    cv = StratifiedKFold(n_splits=_safe_folds(y_train, desired=N_FOLDS),
                         shuffle=True, random_state=RANDOM_STATE)

    cv_acc = cross_val_score(pipe, X_train, y_train, cv=cv, scoring="accuracy")
    # Fit + test
    pipe.fit(X_train, y_train)
    y_pred = pipe.predict(X_test)
    acc = accuracy_score(y_test, y_pred)
    f1 = f1_score(y_test, y_pred)
    # ROC AUC (needs proba or decision)
    try:
        y_proba = pipe.predict_proba(X_test)[:, 1]
    except Exception:
        try:
            y_proba = pipe.decision_function(X_test)
        except Exception:
            y_proba = None
    auc = roc_auc_score(y_test, y_proba) if y_proba is not None else np.nan
    cm = confusion_matrix(y_test, y_pred)

    # Save metrics
    metrics = pd.Series({
        "cv_accuracy_mean": cv_acc.mean(),
        "cv_accuracy_std": cv_acc.std(),
        "test_accuracy": acc,
        "test_f1": f1,
        "test_auc": auc
    })
    metrics.to_csv(outdir / f"{name}__{model_name}__metrics.csv")

    # Confusion matrix plot
    plt.figure()
    im = plt.imshow(cm, cmap="Blues")
    plt.colorbar(im)
    plt.title(f"{name} | {model_name} | Confusion matrix")
    plt.xlabel("Predicted")
    plt.ylabel("True")
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            plt.text(j, i, cm[i, j], ha='center', va='center')
    plt.tight_layout()
    plt.savefig(outdir / f"{name}__{model_name}__confusion_matrix.png", dpi=200)
    plt.close()

    # ROC curve
    if y_proba is not None:
        plt.figure()
        RocCurveDisplay.from_predictions(y_test, y_proba)
        plt.title(f"{name} | {model_name} | ROC")
        plt.tight_layout()
        plt.savefig(outdir / f"{name}__{model_name}__roc.png", dpi=200)
        plt.close()

    return metrics

def run():
    plots = RESULTS / "plots"
    tables = RESULTS / "tables"
    for p in [plots, tables]:
        p.mkdir(parents=True, exist_ok=True)

    X, y, clinical = load_luad_data()

    # ----- Feature set A: Top-K by ANOVA -----
    Xk, feats_k, scores, pvals = _select_topk(X, y)
    Xtr, Xte, ytr, yte = _split(Xk, y)
    bundle_topk = RunBundle("topk", Xtr, Xte, ytr, yte, feats_k)

    # Save feature scores
    pd.DataFrame({
        "feature": feats_k,
        "anova_score": [scores[list(X.columns).index(f)] for f in feats_k],
        "pvalue": [pvals[list(X.columns).index(f)] for f in feats_k]
    }).sort_values("anova_score", ascending=False).to_csv(tables / "topk_feature_scores.csv", index=False)

    # ----- Feature set B: Bio-informed intersection -----
    mirna_prior = read_miRNA_lists()
    Xbio_full, feats_bio = _select_bio(X, mirna_prior)

    if len(feats_bio) == 0:
        print("[INFO] Bio-informed intersection is empty; falling back to all features.")

    # If too many, shrink with ANOVA within the intersection
    if len(feats_bio) > 0 and len(feats_bio) > TOP_K:
        Xbio_full, feats_bio, _, _ = _select_topk(Xbio_full, y, k=TOP_K)

    Xtr_b, Xte_b, ytr_b, yte_b = _split(Xbio_full, y)
    bundle_bio = RunBundle("bio", Xtr_b, Xte_b, ytr_b, yte_b, feats_bio)

    # ----- Train/evaluate -----
    summaries = []
    for bundle in [bundle_topk, bundle_bio]:
        pipes = _build_pipelines(len(bundle.selected_features))
        for mname, pipe in pipes.items():
            metrics = evaluate_and_plot(bundle.name, mname, pipe,
                                        bundle.X_train, bundle.y_train, bundle.X_test, bundle.y_test,
                                        plots)
            row = metrics.to_dict()
            row.update({"feature_set": bundle.name, "model": mname, "n_features": len(bundle.selected_features)})
            summaries.append(row)

    pd.DataFrame(summaries).to_csv(tables / "model_comparison.csv", index=False)

if __name__ == "__main__":
    run()

