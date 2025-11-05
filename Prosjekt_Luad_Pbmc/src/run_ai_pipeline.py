# Orchestrates training and figure generation
from train_models import run as run_models
from eval_plots import heatmap_topk, pca_scatter, rf_feature_importance

if __name__ == "__main__":
    run_models()
    heatmap_topk()
    pca_scatter()
    rf_feature_importance()
    print("Done. See results/ai/plots and results/ai/tables.")
