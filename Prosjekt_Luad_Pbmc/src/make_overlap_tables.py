from pathlib import Path

def read_gene_list(p):
    return [g.strip() for g in Path(p).read_text().splitlines() if g.strip()]

root = Path(__file__).resolve().parents[1]
fg = root / "filtered_gene_lists"
out = root / "results" / "enrichment" / "comparison"
out.mkdir(parents=True, exist_ok=True)

pbmc_up = set(read_gene_list(fg / "pbmc_up.txt"))
pbmc_down = set(read_gene_list(fg / "pbmc_down.txt"))
luad_up = set(read_gene_list(fg / "luad_up.txt"))
luad_down = set(read_gene_list(fg / "luad_down.txt"))

pairs = {
    "exercise_up__luad_down.txt": sorted(pbmc_up & luad_down),
    "exercise_down__luad_up.txt": sorted(pbmc_down & luad_up),
    "exercise_up__luad_up.txt": sorted(pbmc_up & luad_up),
    "exercise_down__luad_down.txt": sorted(pbmc_down & luad_down),
}

for name, genes in pairs.items():
    (out / name).write_text("\n".join(genes))
print("done:", {k: len(v) for k,v in pairs.items()})
