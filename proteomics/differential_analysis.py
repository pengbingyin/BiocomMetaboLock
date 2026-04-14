import json
from pathlib import Path
import numpy as np
from scipy.stats import ttest_ind
from collections import defaultdict
from tqdm import tqdm

print("🚀 Loading gene database for differential analysis...")
DATA_FILE = Path("output") / "gene_database.json"
OUTPUT_FILE = Path("output") / "differential_analysis_results.json"

with open(DATA_FILE, "r", encoding="utf-8") as f:
    db = json.load(f)

print(f"✅ Loaded {len(db)} genes successfully.\n")

# Helper: get valid replicates (ignore None)
def get_replicates(data_dict: dict, condition: str):
    if condition not in data_dict:
        return []
    vals = [v for v in data_dict[condition] if v is not None]
    return vals

# ==================== DEFINE ALL COMPARISONS ====================
session1_comparisons = [  # Within-strain time/condition changes
    {"name": "CENPK_ETH_vs_CENPK_EXP", "condA": "CENPK_ETH", "condB": "CENPK_EXP"},
    {"name": "HUMgE_ETH_vs_HUMgE_EXP", "condA": "HUMgE_ETH", "condB": "HUMgE_EXP"},
    {"name": "HUMgEhE_ETH_vs_HUMgEhE_EXP", "condA": "HUMgEhE_ETH", "condB": "HUMgEhE_EXP"},
    {"name": "HUMgE_72h_vs_HUMgE_ETH", "condA": "HUMgE_72h", "condB": "HUMgE_ETH"},
    {"name": "HUMgEhE_72h_vs_HUMgEhE_ETH", "condA": "HUMgEhE_72h", "condB": "HUMgEhE_ETH"},
    {"name": "NER_96h_vs_NER_72h", "condA": "NER_96h", "condB": "NER_72h"},
    {"name": "NERg_96h_vs_NERg_72h", "condA": "NERg_96h", "condB": "NERg_72h"},
    {"name": "NERgE_96h_vs_NERgE_72h", "condA": "NERgE_96h", "condB": "NERgE_72h"},
    {"name": "VAL_96h_vs_VAL_72h", "condA": "VAL_96h", "condB": "VAL_72h"},
    {"name": "VALg_96h_vs_VALg_72h", "condA": "VALg_96h", "condB": "VALg_72h"},
    {"name": "LIMg_96h_vs_LIMg_72h", "condA": "LIMg_96h", "condB": "LIMg_72h"},
    {"name": "HUMgE_96h_vs_HUMgE_72h", "condA": "HUMgE_96h", "condB": "HUMgE_72h"},
    {"name": "HUMgEhE_96h_vs_HUMgEhE_72h", "condA": "HUMgEhE_96h", "condB": "HUMgEhE_72h"},
]

session2_comparisons = [  # Between-strain at same condition (exclude CENPK)
    {"name": "NERg_72h_vs_NER_72h", "condA": "NERg_72h", "condB": "NER_72h"},
    {"name": "NERgE_72h_vs_NERg_72h", "condA": "NERgE_72h", "condB": "NERg_72h"},
    {"name": "VAL_72h_vs_NER_72h", "condA": "VAL_72h", "condB": "NER_72h"},
    {"name": "VALg_72h_vs_VAL_72h", "condA": "VALg_72h", "condB": "VAL_72h"},
    {"name": "LIMg_72h_vs_NER_72h", "condA": "LIMg_72h", "condB": "NER_72h"},
    {"name": "NERg_96h_vs_NER_96h", "condA": "NERg_96h", "condB": "NER_96h"},
    {"name": "NERgE_96h_vs_NERg_96h", "condA": "NERgE_96h", "condB": "NERg_96h"},
    {"name": "VAL_96h_vs_NER_96h", "condA": "VAL_96h", "condB": "NER_96h"},
    {"name": "VALg_96h_vs_VAL_96h", "condA": "VALg_96h", "condB": "VAL_96h"},
    {"name": "LIMg_96h_vs_NER_96h", "condA": "LIMg_96h", "condB": "NER_96h"},
    {"name": "HUMgE_EXP_vs_HUMg_EXP", "condA": "HUMgE_EXP", "condB": "HUMg_EXP"},
    {"name": "HUMgEh_EXP_vs_HUMgE_EXP", "condA": "HUMgEh_EXP", "condB": "HUMgE_EXP"},
    {"name": "HUMgEhE_EXP_vs_HUMgEh_EXP", "condA": "HUMgEhE_EXP", "condB": "HUMgEh_EXP"},
    {"name": "HUMgEhE_EXP_vs_HUMgE_EXP", "condA": "HUMgEhE_EXP", "condB": "HUMgE_EXP"},


]

session3_comparisons = []
non_cenpk_strains = ["NER", "NERg", "NERgE", "VAL", "VALg", "HUMg", "HUMgE", "HUMgEh", "HUMgEhE", "LIMg"]
for strain in non_cenpk_strains:
    for time in ["72h", "96h", "EXP", "ETH"]:
        cond = f"{strain}_{time}"
        session3_comparisons.append({"name": f"{cond}_vs_CENPK_EXP", "condA": cond, "condB": "CENPK_EXP"})
        if time in ["72h", "96h", "ETH"]:
            session3_comparisons.append({"name": f"{cond}_vs_CENPK_ETH", "condA": cond, "condB": "CENPK_ETH"})

all_comparisons = session1_comparisons + session2_comparisons + session3_comparisons

print(f"Session 1 (within-strain): {len(session1_comparisons)} comparisons")
print(f"Session 2 (between-strain): {len(session2_comparisons)} comparisons")
print(f"Session 3 (vs CENPK): {len(session3_comparisons)} comparisons")
print(f"Total comparisons to run: {len(all_comparisons)}\n")

# ====================== RUN ANALYSIS WITH DETAILED PROGRESS ======================
results = defaultdict(dict)

print("🔄 Starting differential analysis...\n")

gene_list = list(db.keys())
for i, gene in enumerate(tqdm(gene_list, desc="Processing genes", unit="gene")):
    linear_data = db[gene].get("linear_data", {})
    log2_data = db[gene].get("log2_data", {})

    gene_comparisons_count = 0

    for comp in all_comparisons:
        name = comp["name"]
        condA = comp["condA"]
        condB = comp["condB"]

        lin_A = [v for v in linear_data.get(condA, []) if v is not None]
        lin_B = [v for v in linear_data.get(condB, []) if v is not None]
        log_A = [v for v in log2_data.get(condA, []) if v is not None]
        log_B = [v for v in log2_data.get(condB, []) if v is not None]

        if len(lin_A) < 2 or len(lin_B) < 2 or len(log_A) < 2 or len(log_B) < 2:
            continue

        # Fold change from log2 means
        mean_log_A = np.mean(log_A)
        mean_log_B = np.mean(log_B)
        foldchange = mean_log_A - mean_log_B

        # Welch's t-test (two-tailed)
        _, p_linear = ttest_ind(lin_A, lin_B, equal_var=False, nan_policy='omit')
        _, p_log2   = ttest_ind(log_A, log_B, equal_var=False, nan_policy='omit')

        results[gene][name] = {
            "foldchange": round(float(foldchange), 4),
            "p_value_linear": round(float(p_linear), 6),
            "p_value_log2": round(float(p_log2), 6)
        }
        gene_comparisons_count += 1

    # Detailed progress print every 200 genes + first 5 genes
    if i < 5 or (i + 1) % 200 == 0:
        print(f"  → Processed gene {i+1:4d}/{len(gene_list)}: {gene} | {gene_comparisons_count} valid comparisons")

print("\n✅ All genes processed successfully!")

# ====================== SAVE RESULTS ======================
with open(OUTPUT_FILE, 'w', encoding='utf-8') as f:
    json.dump(results, f, indent=2)

print(f"\n📁 Results saved to: {OUTPUT_FILE.resolve()}")
print(f"   📊 Genes with at least one comparison: {len(results)}")
print(f"   📋 Total comparison entries: {sum(len(v) for v in results.values())}")
print("\n🎉 Analysis complete! Open output/differential_analysis_results.json to view results.")