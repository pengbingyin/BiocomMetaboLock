import pandas as pd
import os

# ====================== PATHS ======================
INPUT_CSV = "data/peptides.csv"
OUTPUT_DIR = "output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

OUTPUT_EXCEL = os.path.join(OUTPUT_DIR, "protein_quantification_proteotypic.xlsx")
OUTPUT_CSV   = os.path.join(OUTPUT_DIR, "protein_quantification_proteotypic.csv")

# ====================== LOAD DATA ======================
print("Loading peptide CSV...")
df = pd.read_csv(INPUT_CSV, low_memory=False)

print(f"Original peptides: {len(df):,}")

# === FIXED FILTER: Use string "TRUE" instead of boolean True ===
df_proteo = df[df['PEP.IsProteotypic'].astype(str).str.upper() == "TRUE"].copy()

print(f"Proteotypic peptides kept: {len(df_proteo):,}")

if len(df_proteo) == 0:
    print("⚠️  No proteotypic peptides found! Check the exact values in column 'PEP.IsProteotypic'")
    print(df['PEP.IsProteotypic'].value_counts(dropna=False))
    exit()

# Identify quantity columns
qty_cols = [col for col in df.columns if "EG.TotalQuantity" in col]

# Clean data: replace "Filtered" and NaN with 0
print("Cleaning quantity data...")
for col in qty_cols:
    df_proteo[col] = df_proteo[col].replace(['Filtered', 'NaN', '', 'nan'], 0)
    df_proteo[col] = pd.to_numeric(df_proteo[col], errors='coerce').fillna(0)

# ====================== SUM TO PROTEIN ======================
print("Summing proteotypic peptides to protein groups...")

protein_df = df_proteo.groupby('PG.ProteinAccessions', as_index=False).agg({
    'PG.ProteinDescriptions': 'first',
    'PG.FastaFiles': 'first',
    **{col: 'sum' for col in qty_cols}
})

protein_df['Num_Peptides'] = df_proteo.groupby('PG.ProteinAccessions').size().values

print(f"Final protein table: {protein_df.shape[0]:,} proteins")

# ====================== SAVE ======================
protein_df.to_excel(OUTPUT_EXCEL, index=False)
protein_df.to_csv(OUTPUT_CSV, index=False)

print("\n✅ Done! Proteotypic-only protein quantification saved.")
print(f"   Excel: {OUTPUT_EXCEL}")
print(f"   CSV:   {OUTPUT_CSV}")