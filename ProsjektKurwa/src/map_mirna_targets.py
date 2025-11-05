# src/map_mirna_targets.py
import os
from collections import defaultdict
import gseapy as gp
import pandas as pd

BASE = os.path.dirname(__file__)
IN_DIR = os.path.join(BASE, "..", "filtered_gene_lists")
OUT_DIR = os.path.join(BASE, "..", "results", "enrichment", "comparison")
os.makedirs(OUT_DIR, exist_ok=True)

luad_up_file = os.path.join(IN_DIR, "luad_up.txt")
luad_down_file = os.path.join(IN_DIR, "luad_down.txt")

def load_list(path):
    with open(path, "r", encoding="utf-8") as f:
        return [ln.strip() for ln in f if ln.strip()]

luad_up_mirs = load_list(luad_up_file)
luad_down_mirs = load_list(luad_down_file)

# ----- fetch miRTarBase gene-set library (term -> list of gene symbols) -----
print("downloading miRTarBase_2017 library ...")
lib = gp.get_library(name="miRTarBase_2017", organism="Human")  # dict: term -> [GENE,...]
lib_terms = set(lib.keys())

# optional: peek at a few keys to see the style
peek_terms = list(lib_terms)[:20]
peek_path = os.path.join(OUT_DIR, "DEBUG_miRTarBase_terms_preview.txt")
with open(peek_path, "w") as f:
    f.write("\n".join(peek_terms))
print(f"saved a preview of library terms to: {peek_path}")

# ----- normalization helpers -----
def strip_hsa(x: str) -> str:
    return x[4:] if x.lower().startswith("hsa-") else x

def strip_arm(x: str) -> str:
    # remove '-3p' or '-5p' suffixes if present
    for suf in ("-3p", "-5p", "_3p", "_5p"):
        if x.lower().endswith(suf):
            return x[: -len(suf)]
    return x

def canon_mir(x: str) -> str:
    # unify case and dash style for robust comparison
    x = x.strip()
    x = x.replace("MIR", "miR").replace("Mir", "miR")
    x = x.replace("_", "-")
    return x

def variants(miR: str):
    """generate possible variants that may exist in the library"""
    a = canon_mir(miR)
    b = canon_mir(strip_hsa(a))
    c = canon_mir(strip_arm(b))
    # also try removing both arm and hsa in different orders
    d = canon_mir(strip_arm(a))
    e = canon_mir(strip_hsa(d))
    # deduplicate while preserving order
    seen, out = set(), []
    for v in [a, b, c, d, e]:
        if v not in seen:
            out.append(v); seen.add(v)
    return out

# build a reverse index: for quick lookup by canonical forms
index = defaultdict(list)
for term in lib_terms:
    t = canon_mir(term)
    index[t].append(term)
    # also index arm-less (if present)
    ta = canon_mir(strip_arm(t))
    index[ta].append(term)

def map_one(miR: str):
    """return list of genes for the first matching library term, with a priority on exact, then arm-less."""
    for v in variants(miR):
        # priority 1: exact canonical match
        if v in lib:
            return lib[v], v, v
        # priority 2: via index (maps to actual library keys)
        if v in index:
            # choose a preferred key: exact first, otherwise first in list
            keys = index[v]
            # try to pick an exact-arm key if exists
            preferred = None
            for k in keys:
                if canon_mir(k) == v:
                    preferred = k; break
            if preferred is None:
                preferred = keys[0]
            return lib[preferred], v, preferred
    return [], None, None

def map_list(mirs, label):
    targets = set()
    unmatched = []
    matches = []
    for m in mirs:
        genes, query_used, lib_key = map_one(m)
        if genes:
            targets.update(g.strip() for g in genes if g and g.strip() and g.strip().upper() != "NA")
            matches.append((m, query_used, lib_key, len(genes)))
        else:
            unmatched.append(m)
    # save logs
    log_ok = os.path.join(OUT_DIR, f"{label}_miRNA_match_log.csv")
    if matches:
        pd.DataFrame(matches, columns=["input_miR", "query_variant", "library_key", "n_genes"]).to_csv(log_ok, index=False)
    log_bad = os.path.join(OUT_DIR, f"{label}_miRNA_unmatched.txt")
    if unmatched:
        with open(log_bad, "w") as f:
            f.write("\n".join(unmatched))
    return sorted(targets), log_ok, log_bad

up_targets, log_up_ok, log_up_bad = map_list(luad_up_mirs, "LUAD_up")
down_targets, log_down_ok, log_down_bad = map_list(luad_down_mirs, "LUAD_down")

# write targets
up_out = os.path.join(OUT_DIR, "luad_up_targets.txt")
down_out = os.path.join(OUT_DIR, "luad_down_targets.txt")
with open(up_out, "w") as f:
    f.write("\n".join(up_targets))
with open(down_out, "w") as f:
    f.write("\n".join(down_targets))

print(f"LUAD-up miRNAs: {len(luad_up_mirs)}  -> mapped targets: {len(up_targets)}")
print(f"LUAD-down miRNAs: {len(luad_down_mirs)} -> mapped targets: {len(down_targets)}")
print("saved target lists:")
print(up_out)
print(down_out)
print("match logs:")
print(log_up_ok, "(if exists)")
print(log_down_ok, "(if exists)")
print("unmatched lists (if any):")
print(log_up_bad)
print(log_down_bad)
