import pandas as pd
import json
from pathlib import Path
from collections import defaultdict

# ====================== CONFIG ======================
DATA_DIR = Path("data")
OUTPUT_DIR = Path("output")
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

NORMALISED_FILE = DATA_DIR / "proteomics_normalised.xlsx"
LOG_FILE = DATA_DIR / "proteomics_log.xlsx"
REGULATORY_FILE = DATA_DIR / "regulatory matrix.xlsx"

OUTPUT_JSON = OUTPUT_DIR / "gene_database.json"

print("🚀 Starting gene database creation...")

# Load proteomics data (main source - we keep ALL genes)
norm_df = pd.read_excel(NORMALISED_FILE, sheet_name="Sheet1", index_col=0)
log_df = pd.read_excel(LOG_FILE, sheet_name="Sheet1", index_col=0)

# Load regulatory matrix
reg_df = pd.read_excel(REGULATORY_FILE, sheet_name="Sheet1", index_col=0)
reg_df = reg_df.fillna(0).astype(int)

# Clean indices
norm_df.index = norm_df.index.astype(str).str.strip()
log_df.index = log_df.index.astype(str).str.strip()
reg_df.index = reg_df.index.astype(str).str.strip()

print(f"Loaded {len(norm_df)} genes from proteomics files")
print(f"Loaded {len(reg_df)} genes from regulatory matrix")

all_genes = sorted(norm_df.index)
print(f"Processing ALL {len(all_genes)} genes from proteomics...")

# Get all transcription factors
all_trans_factors = list(reg_df.columns) if not reg_df.empty else []

database = {}

for gene in all_genes:
    linear_dict = defaultdict(lambda: [None, None, None])
    log2_dict = defaultdict(lambda: [None, None, None])

    # Parse sample columns: Strain_Replicate_Time → Condition = Strain_Time
    for col in norm_df.columns:
        parts = str(col).split('_')
        if len(parts) >= 3:
            strain = parts[0]
            try:
                rep_idx = int(parts[1]) - 1
                time_part = '_'.join(parts[2:])
                condition = f"{strain}_{time_part}"   # e.g. "CENPK_EXP", "NER_72h"

                if 0 <= rep_idx <= 2:
                    linear_dict[condition][rep_idx] = float(norm_df.loc[gene, col])
                    log2_dict[condition][rep_idx]   = float(log_df.loc[gene, col])
            except (ValueError, IndexError, KeyError):
                continue

    # === NEW RULE for regulation ===
    if ';' in gene:
        # Protein group (e.g. TEF1;TEF2) → set regulation to empty dict
        regulation = {}
        group_note = " (protein group → regulation set to empty)"
    else:
        # Single gene → take regulation from matrix or empty if not present
        if gene in reg_df.index:
            regulation = reg_df.loc[gene].to_dict()
        else:
            regulation = {}
        group_note = ""

    database[gene] = {
        "linear_data": dict(linear_dict),
        "log2_data":   dict(log2_dict),
        "regulation":  regulation
    }

# ====================== SAVE TO OUTPUT FOLDER ======================
with open(OUTPUT_JSON, 'w', encoding='utf-8') as f:
    json.dump(database, f, indent=2, ensure_ascii=False)

print(f"\n✅ Success! Gene database saved.")
print(f"   📁 Location        : {OUTPUT_JSON.resolve()}")
print(f"   📊 Total genes     : {len(database)}")
print(f"   📋 Protein groups (with empty regulation): {sum(1 for g in database if ';' in g)}")

# Example output
if database:
    example_gene = next(iter(database))
    ex = database[example_gene]
    is_group = ';' in example_gene
    print(f"\nExample gene: '{example_gene}'{ ' (protein group → regulation empty)' if is_group else '' }")
    print(f"   Conditions : {len(ex['linear_data'])}")
    print(f"   Regulators : {len(ex['regulation'])} {'(empty for protein groups)' if is_group else ''}")
    print(f"   First few conditions: {list(ex['linear_data'].keys())[:5]} ...")