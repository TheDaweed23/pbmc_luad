from pathlib import Path

def clean_lines(p: Path):
    """
    Clean up a text file containing gene names.
    Steps:
    1. Read all lines from the file.
    2. Strip whitespace and surrounding quotes.
    3. Ignore empty lines.
    4. If a line contains multiple gene names separated by "///",
       keep only the first one.
    5. Remove duplicates and sort the final list alphabetically.
    6. Overwrite the file with the cleaned list.
    """
    out = []

    # Read and process each line
    for ln in p.read_text(encoding="utf-8").splitlines():
        ln = ln.strip().strip('"')  # remove leading/trailing spaces and quotes
        if not ln:
            continue  # skip empty lines

        # Keep only the first token if multiple gene names appear
        ln = ln.split("///")[0].strip()
        out.append(ln)

    # Remove duplicates, sort, and write back to the same file
    p.write_text("\n".join(sorted(set(out))), encoding="utf-8")


# Define project root (parent of the current script’s directory)
root = Path(__file__).resolve().parents[1]

# Clean up both PBMC gene lists
for name in ["pbmc_up.txt", "pbmc_down.txt"]:
    clean_lines(root / "filtered_gene_lists" / name)

print("PBMC lists cleaned.")
