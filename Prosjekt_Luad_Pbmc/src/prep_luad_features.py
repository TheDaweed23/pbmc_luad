# src/prep_luad_features.py

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Tuple, List, Dict, Optional

from ai_config import LUAD_XLSX, LUAD_UP_TXT, LUAD_DOWN_TXT, RESULTS


def _read_excel_sheets(xlsx_path: Path) -> Dict[str, pd.DataFrame]:
    """
    Open the LUAD Excel and return all sheets as {sheet_name: DataFrame}.
    """
    try:
        xl = pd.ExcelFile(xlsx_path)
    except FileNotFoundError as e:
        raise FileNotFoundError(
            f"Excel not found at: {xlsx_path}\n"
            f"Tip: place 'LUAD_expression_with_clinical.xlsx' in data/raw/ or data/"
        ) from e
    return {name: xl.parse(name) for name in xl.sheet_names}


def _detect_sheets(sheets: Dict[str, pd.DataFrame]) -> Tuple[pd.DataFrame, Optional[pd.DataFrame]]:
    """
    Try to find the expression sheet (has 'miRBaseName') and a clinical-like sheet
    (has columns such as 'Sample', 'Age', 'Gender').
    Returns: (expr_df, clinical_df or None)
    """
    expr = None
    clinical = None
    for name, df in sheets.items():
        cols_lower = [str(c).lower() for c in df.columns]
        if "mirbasename" in cols_lower:
            expr = df.copy()
        if {"sample", "age", "gender"}.issubset(set(cols_lower)):
            clinical = df.copy()
    if expr is None:
        raise ValueError("Could not find an expression sheet with a 'miRBaseName' column.")
    return expr, clinical


def load_luad_data() -> Tuple[pd.DataFrame, pd.Series, pd.DataFrame]:
    """
    Load miRNA expression and labels from the LUAD Excel.

    Returns:
      X: samples x features DataFrame (numeric), contains both LUAD_* and HC_* samples
      y: Series of labels (1=LUAD, 0=HC) indexed by sample name
      clinical: clinical DataFrame (may contain NaNs for HC_* rows)
    """
    sheets = _read_excel_sheets(LUAD_XLSX)
    expr, clinical = _detect_sheets(sheets)

    # Normalize the feature-name column, set as index
    expr = expr.rename(columns={expr.columns[0]: "miRBaseName"})
    expr.set_index("miRBaseName", inplace=True)

    # Build labels from column prefixes
    cols = expr.columns.astype(str)
    is_luad = cols.str.startswith("LUAD_")
    is_hc = cols.str.startswith("HC_")
    if not (is_luad.any() and is_hc.any()):
        raise ValueError(
            "Expected expression columns starting with 'LUAD_' and 'HC_'. "
            f"Found columns like: {list(cols[:6])} ..."
        )
    y = pd.Series(np.where(is_luad, 1, 0), index=cols, name="label")

    # Transpose to samples x features and coerce to numeric
    X = expr.T
    X = X.apply(pd.to_numeric, errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(0.0)

    # Drop zero-variance features (often all-zero miRNAs)
    var = X.var(axis=0)
    keep = var[var > 0].index
    if len(keep) == 0:
        raise ValueError("All features have zero variance after cleaning.")
    X = X[keep]

    # Attach clinical without filtering out HC samples (HC will be NaN in clinical)
    if clinical is not None and "Sample" in clinical.columns:
        clinical = clinical.set_index("Sample")
        clinical = clinical.reindex(X.index)
    else:
        clinical = pd.DataFrame(index=X.index)

    # Save quick summary
    (RESULTS / "artifacts").mkdir(parents=True, exist_ok=True)
    X.describe().T.to_csv(RESULTS / "artifacts" / "miRNA_feature_summary.csv")

    return X, y, clinical


def read_miRNA_lists() -> List[str]:
    """
    Read LUAD up/down lists (miRNA IDs). Returns a de-duplicated list of IDs.
    If files are missing, returns an empty list (the calling code handles fallback).
    """
    ids: List[str] = []
    for path in [LUAD_UP_TXT, LUAD_DOWN_TXT]:
        try:
            with open(path, "r", encoding="utf-8") as f:
                for line in f:
                    t = line.strip()
                    if t:
                        ids.append(t)
        except FileNotFoundError:
            # ignore; caller will fallback if intersection becomes empty
            pass

    # De-duplicate while preserving order
    seen = set()
    out: List[str] = []
    for x in ids:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out
