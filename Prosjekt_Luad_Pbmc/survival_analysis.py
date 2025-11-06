import pandas as pd
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import matplotlib.pyplot as plt
import os

# 1. Load the cleaned data
df = pd.read_excel("LUAD_cleaned.xlsx")

# 2. Define the column names for time and event status in the survival analysis.
time_col = "OS_days"  # Overall Survival time in days for each patient
status_col = None  # We do not have an explicit event/censor status column (assuming all events or using OS_days as time)

# Identify all miRNA columns.
# Based on how we constructed LUAD_cleaned.xlsx, clinical columns are first (Sample, Age, etc.), and miRNAs start from column 6 onward.
miRNA_cols = df.columns[
    6:]  # columns after 'Smoking_Status' (adjust index if needed based on actual columns in cleaned data)

# 3. Create an output directory for survival plots
os.makedirs("Survival_Plots", exist_ok=True)

# 4. Prepare a list to collect results (miRNA names and p-values)
results = []

# 5. Loop through each miRNA and perform Kaplan-Meier analysis
for mir in miRNA_cols:
    try:
        # Split patients into high and low expression groups by median value of this miRNA
        median_value = df[mir].median()
        high_group = df[df[mir] > median_value]
        low_group = df[df[mir] <= median_value]

        # Ensure groups are not empty (in case all values are identical, skip)
        if len(high_group) == 0 or len(low_group) == 0:
            continue

        # Initialize Kaplan-Meier fitters for each group
        kmf_high = KaplanMeierFitter()
        kmf_low = KaplanMeierFitter()

        # Fit the Kaplan-Meier model on survival time for each group.
        # event_observed is set to None, meaning we assume all subjects reached the event (no censoring data available).
        kmf_high.fit(durations=high_group[time_col], event_observed=None, label=f"High {mir}")
        kmf_low.fit(durations=low_group[time_col], event_observed=None, label=f"Low {mir}")

        # Perform the log-rank test to compare survival between high vs low groups
        # Since we don't have censoring info, we provide the durations and assume all events.
        lr_result = logrank_test(high_group[time_col], low_group[time_col],
                                 event_observed_A=None, event_observed_B=None)
        p_value = lr_result.p_value

        # Record the result
        results.append({"miRNA": mir, "p_value": p_value})

        # 6. If the difference is statistically significant (p < 0.05), save a KM plot for this miRNA
        if p_value < 0.05:
            plt.figure()  # create new figure
            # Plot the KM curves for both groups on the same axes
            kmf_high.plot(ci_show=False)  # ci_show=False to hide confidence intervals for clarity
            kmf_low.plot(ci_show=False)
            # Add title and labels
            plt.title(f"{mir} (Log-rank p = {p_value:.4f})")
            plt.xlabel("Days")
            plt.ylabel("Survival Probability")
            plt.tight_layout()
            # Save plot as PNG in the Survival_Plots directory
            plt.savefig(f"Survival_Plots/{mir}.png")
            plt.close()  # close figure to avoid overlap in next plots

    except Exception as e:
        # If any error occurs (for example, if data is insufficient), skip this miRNA
        print(f"Skipped {mir} due to error: {e}")
        continue

# 7. After loop, compile results into a DataFrame and save the summary
res_df = pd.DataFrame(results).sort_values("p_value")
res_df.to_excel("Significant_miRNAs.xlsx", index=False)  # save all miRNAs with their p-values to Excel
res_df.to_csv("Significant_miRNAs.csv", index=False)  # also save as CSV for convenience

# Also print top results to console for quick viewing
print("\nTop 10 miRNAs by significance (lowest p-values):")
print(res_df.head(10))

print(f"\n Survival analysis complete! {len(results)} miRNAs analyzed.")
print(
    "Results saved in 'Significant_miRNAs.xlsx/CSV'. Kaplan-Meier plots for p<0.05 are in the 'Survival_Plots' folder.")
