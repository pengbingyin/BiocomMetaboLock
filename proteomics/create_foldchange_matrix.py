import json
from pathlib import Path
import pandas as pd
import numpy as np

print("🚀 Loading differential analysis results...")

# Input and output paths
INPUT_JSON = Path("output") / "differential_analysis_results.json"
OUTPUT_CSV = Path("output") / "foldchange_matrix.csv"

# Load the JSON
with open(INPUT_JSON, "r", encoding="utf-8") as f:
    data = json.load(f)

print(f"✅ Loaded results for {len(data)} genes.")

# Build the fold change matrix
print("🔄 Building fold change matrix...")

# Collect all unique comparison names
all_comparisons = set()
for gene_comps in data.values():
    all_comparisons.update(gene_comps.keys())

comparison_list = sorted(all_comparisons)   # sort for consistent columns

print(f"   Found {len(comparison_list)} unique comparisons.")

# Create matrix: gene → comparison → foldchange
matrix = {}
for gene, comps in data.items():
    row = {}
    for comp in comparison_list:
        # Get foldchange if the comparison exists for this gene, otherwise NaN
        foldchange = comps.get(comp, {}).get("foldchange")
        row[comp] = foldchange if foldchange is not None else np.nan
    matrix[gene] = row

# Convert to pandas DataFrame
df = pd.DataFrame.from_dict(matrix, orient='index')
df = df[comparison_list]        # enforce column order
df = df.round(4)                # round to 4 decimal places

print(f"✅ Matrix created successfully!")
print(f"   Shape: {df.shape[0]} genes × {df.shape[1]} comparisons")

# Save to CSV
df.to_csv(OUTPUT_CSV, index=True)   # index=True keeps gene names as first column

print(f"\n📁 Fold change matrix saved to:")
print(f"   {OUTPUT_CSV.resolve()}")
print(f"   File size: {OUTPUT_CSV.stat().st_size / (1024*1024):.2f} MB")

# Preview
print("\n📊 Preview (first 5 genes × first 8 comparisons):")
preview = df.iloc[:5, :8]
print(preview)

print("\n🎉 Done! You can now open 'foldchange_matrix.csv' in Excel.")