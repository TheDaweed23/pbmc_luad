# Paths and knobs for the AI
from pathlib import Path

# project-root relative
ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results" / "ai"
RESULTS.mkdir(parents=True, exist_ok=True)

# input files
LUAD_XLSX = DATA / "raw" / "LUAD_expression_with_clinical.xlsx"
LUAD_UP_TXT = DATA / "filtered_gene_lists" / "luad_up.txt"     # miRNA IDs
LUAD_DOWN_TXT = DATA / "filtered_gene_lists" / "luad_down.txt" # miRNA IDs

# feature selection
TOP_K = 100                 # number of top features by ANOVA F-test
USE_BIO_OVERLAP = True      # also build a "bio-informed" feature set (intersection with up/down lists)

# models to train (two taught + one exploratory)
MODELS = ["logreg", "knn", "rf"]   # logistic regression, kNN, random forest

# cv/eval
TEST_SIZE = 0.25
RANDOM_STATE = 42
N_FOLDS = 5
