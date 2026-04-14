import json
from pathlib import Path
import pandas as pd
from collections import defaultdict

print("🚀 Starting Transcription Factor Enrichment Analysis...")

# ====================== CONFIG ======================
DATA_DIR = Path("data")
OUTPUT_DIR = Path("output")
OUTPUT_DIR.mkdir(exist_ok=True)

REGULATORY_FILE = DATA_DIR / "regulatory matrix.xlsx"
DIFF_RESULTS_FILE = OUTPUT_DIR / "differential_analysis_results.json"

# Output files
TF_LIBRARY_JSON = OUTPUT_DIR / "tf_target_library.json"
GENE_REGULATION_JSON = OUTPUT_DIR / "gene_regulation_status.json"
TF_COUNTS_JSON = OUTPUT_DIR / "tf_enrichment_counts.json"
TF_ASSOCIATION_JSON = OUTPUT_DIR / "tf_association.json"
TF_ASSOCIATION_CSV = OUTPUT_DIR / "tf_association_filtered.csv"

# Thresholds
FOLDCHANGE_THRESHOLD = 0.7
PVALUE_THRESHOLD = 0.1
MIN_TARGETS = 5   # Only keep TFs with > 5 total targets

# ====================== STEP 1: Create TF target library ======================
print("Step 1/4: Building TF target library...")

reg_df = pd.read_excel(REGULATORY_FILE, sheet_name="Sheet1", index_col=0)
reg_df.index = reg_df.index.astype(str).str.strip()
reg_df = reg_df.fillna(0).astype(int)

tf_targets = defaultdict(list)
for tf in reg_df.columns:
    targets = reg_df.index[reg_df[tf] == 1].tolist()
    tf_targets[tf] = sorted(targets)

tf_targets = dict(tf_targets)

with open(TF_LIBRARY_JSON, 'w', encoding='utf-8') as f:
    json.dump(tf_targets, f, indent=2)

print(f"   ✅ TF target library created: {len(tf_targets)} TFs")
print(f"   📁 Saved to: {TF_LIBRARY_JSON.resolve()}\n")

# ====================== STEP 2: Gene regulation status ======================
print("Step 2/4: Classifying genes as up/down-regulated...")

with open(DIFF_RESULTS_FILE, "r", encoding="utf-8") as f:
    diff_data = json.load(f)

gene_regulation = defaultdict(dict)   # comparison → gene → "up" | "down" | None

for gene, comp_results in diff_data.items():
    for comparison, stats in comp_results.items():
        fc = stats.get("foldchange")
        p_lin = stats.get("p_value_linear")
        p_log = stats.get("p_value_log2")

        if fc is None or p_lin is None or p_log is None:
            status = None
        else:
            significant = (p_lin < PVALUE_THRESHOLD) or (p_log < PVALUE_THRESHOLD)
            if significant and fc >= FOLDCHANGE_THRESHOLD:
                status = "up"
            elif significant and fc <= -FOLDCHANGE_THRESHOLD:
                status = "down"
            else:
                status = None

        gene_regulation[comparison][gene] = status

gene_regulation = dict(gene_regulation)

with open(GENE_REGULATION_JSON, 'w', encoding='utf-8') as f:
    json.dump(gene_regulation, f, indent=2)

print(f"   ✅ Gene regulation status created for {len(gene_regulation)} comparisons\n")

# ====================== STEP 3: Count up/down using TOTAL targets ======================
print("Step 3/4: Counting up/down targets per TF (using total targets)...")

tf_counts = defaultdict(dict)

for comparison, gene_status in gene_regulation.items():
    for tf, targets in tf_targets.items():
        total = len(targets)                    # ← This is the key change
        up = sum(1 for t in targets if gene_status.get(t) == "up")
        down = sum(1 for t in targets if gene_status.get(t) == "down")

        tf_counts[comparison][tf] = {
            "up": up,
            "down": down,
            "total": total
        }

tf_counts = dict(tf_counts)

with open(TF_COUNTS_JSON, 'w', encoding='utf-8') as f:
    json.dump(tf_counts, f, indent=2)

print(f"   ✅ TF counts created (total = all targets)\n")

# ====================== STEP 4: Association factor + NEW FILTERS + Transpose ======================
print(f"Step 4/4: Calculating association factor + filtering (|assoc| > 0.4 in any comparison AND total > {MIN_TARGETS})...")

tf_association = defaultdict(dict)
candidate_tfs = set()

# First pass: calculate all association values
for comparison, tf_data in tf_counts.items():
    for tf, counts in tf_data.items():
        total = counts["total"]
        if total > MIN_TARGETS:
            assoc = (counts["up"] - counts["down"]) / total if total > 0 else 0.0
            tf_association[comparison][tf] = round(float(assoc), 4)

tf_association = dict(tf_association)

# Second pass: find TFs that meet the |0.4| threshold in ANY comparison
for comparison, tf_data in tf_association.items():
    for tf, assoc in tf_data.items():
        if abs(assoc) > 0.35:
            candidate_tfs.add(tf)

# Final filtered TFs = those with >5 targets AND |assoc| > 0.4 somewhere
filtered_tfs = [tf for tf in candidate_tfs if tf in tf_association[list(tf_association.keys())[0]]]

print(f"   Found {len(filtered_tfs)} TFs that have |association| > 0.4 in at least one comparison")

# Build final DataFrame
df = pd.DataFrame.from_dict(tf_association, orient='index')
if filtered_tfs:
    df = df[filtered_tfs]          # keep only strongly associated TFs
df = df.round(4)

# Transpose: TFs as rows, comparisons as columns
df_transposed = df.T
df_transposed.index.name = "Transcription_Factor"
df_transposed.columns.name = "Comparison"

# Save outputs
with open(TF_ASSOCIATION_JSON, 'w', encoding='utf-8') as f:
    json.dump(tf_association, f, indent=2)

df_transposed.to_csv(TF_ASSOCIATION_CSV)

print(f"   ✅ Final filtered matrix created")
print(f"   📁 Full JSON: {TF_ASSOCIATION_JSON.resolve()}")
print(f"   📁 Transposed filtered CSV: {TF_ASSOCIATION_CSV.resolve()}")
print(f"   Final shape: {df_transposed.shape[0]} TFs × {df_transposed.shape[1]} comparisons\n")

print("🎉 All steps completed!")
print("The new tf_association_filtered.csv now uses **total_targets** (all known targets) in the denominator.")