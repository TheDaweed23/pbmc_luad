import pandas as pd

# 1. Load the Excel file with two sheets: "Expression" and "Clinical"
xls = pd.ExcelFile("LUAD_expression_with_clinical.xlsx")  # ensure this file path is correct
expr_df = pd.read_excel(xls, sheet_name="Expression")
clin_df = pd.read_excel(xls, sheet_name="Clinical")

# 2. Filter expression data to keep only LUAD samples.
# We assume sample columns start with "LUAD_" for LUAD patients.
luad_sample_cols = ["miRBaseName"] + [col for col in expr_df.columns if col.startswith("LUAD_")]

# Subset the expression DataFrame to only include LUAD samples and the miRNA name column
luad_expr = expr_df[luad_sample_cols]

# 3. Transpose the expression DataFrame so that each row is a sample and each column is a miRNA.
# Currently, rows are miRNAs and columns are samples; we want rows = samples.
luad_expr = luad_expr.set_index("miRBaseName").T  # set miRNA names as index, then transpose
luad_expr.reset_index(inplace=True)              # make the sample IDs a column again
luad_expr = luad_expr.rename(columns={"index": "Sample"})  # rename the index column to "Sample"

# 4. Merge the transposed expression data with the clinical data on the "Sample" column.
# This will add clinical columns (age, gender, OS_days, etc.) alongside each sample's miRNA data.
merged_data = pd.merge(clin_df, luad_expr, on="Sample")

# 5. Save the cleaned and merged data to a new Excel file for use in analysis.
merged_data.to_excel("LUAD_cleaned.xlsx", index=False)
print(" Data cleaning complete! The merged data is saved as 'LUAD_cleaned.xlsx'.")
