import pandas as pd
import numpy as np
import os
import matplotlib

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# ====================== CONFIG ======================
INPUT_FILE = "data/input_PCA.xlsx"
OUTPUT_DIR = "output/pca_output_percent"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ====================== LOAD & PREPROCESS ======================
df = pd.read_excel(INPUT_FILE, sheet_name="Sheet1")
sample_cols = df.columns[1:]

X = df[sample_cols].T.copy()
X = X.replace(0, np.nan)
X = np.log2(X + 1)
X = X.fillna(0)

print(f"Loaded {df.shape[0]} proteins × {len(sample_cols)} samples")

# ====================== PARSE METADATA ======================
metadata = []
for sample in X.index:
    parts = str(sample).split('_')
    if len(parts) >= 3:
        strain = parts[0]
        replicate = parts[1]
        condition = '_'.join(parts[2:])
    else:
        strain = sample
        replicate = "NA"
        condition = "NA"

    metadata.append({
        "Sample": sample,
        "Strain": strain,
        "Replicate": replicate,
        "Condition": condition
    })

meta_df = pd.DataFrame(metadata)

# ====================== PCA ======================
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
pca = PCA(n_components=10)
principal_components = pca.fit_transform(X_scaled)

pca_df = pd.DataFrame(principal_components, columns=[f"PC{i + 1}" for i in range(10)], index=X.index)
pca_df = pca_df.join(meta_df.set_index("Sample"))

explained = pca.explained_variance_ratio_ * 100
print("\nExplained variance (%):")
for i, v in enumerate(explained[:5], 1):
    print(f"PC{i}: {v:.2f}%")

# ====================== COLORS ======================
condition_colors = {
    'EXP': '#ff7f0e',
    'ETH': '#d62728',
    '72h': '#1f77b4',
    '96h': '#0d3b66'
}

# ====================== 1. MAIN FIGURE ======================
width_in = 17 / 2.54
fig, ax = plt.subplots(figsize=(width_in, width_in))

sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="Condition", style="Strain",
                s=80, alpha=1.0, palette=condition_colors, edgecolor="black", linewidth=0.8, ax=ax)

# Replicate numbers
for _, row in pca_df.iterrows():
    ax.text(row["PC1"], row["PC2"] + 1.2, str(row["Replicate"]),
            fontsize=10, color="grey", ha='center', va='center',
            fontweight='bold', fontfamily='Arial')

ax.set_title("PCA on Percentage-to-Average Normalized Proteomics Data", fontsize=13, pad=20)
ax.set_xlabel(f"PC1 ({explained[0]:.1f}% variance)")
ax.set_ylabel(f"PC2 ({explained[1]:.1f}% variance)")

if ax.get_legend() is not None:
    ax.get_legend().remove()

plt.tight_layout()
fig.savefig(os.path.join(OUTPUT_DIR, "PCA_main_figure_17cm_s80.png"), dpi=300, bbox_inches='tight')

# ====================== 2. GROUPED LEGEND + EXTRACT MAPPING ======================
print("\n=== Extracting marker mapping from grouped legend ===")

legend_fig, legend_ax = plt.subplots(figsize=(14, 4))
dummy = sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="Condition", style="Strain",
                        s=80, alpha=1.0, palette=condition_colors, edgecolor="black", linewidth=0.8, ax=legend_ax)

handles, labels = dummy.get_legend_handles_labels()
legend_ax.clear()
legend_ax.axis('off')

# Print the mapping for debugging
print("Strain -> Marker mapping:")
strain_to_marker = {}
for h, label in zip(handles, labels):
    if "Strain" in label or label == "Strain":
        continue
    # Try to extract strain from label
    if "—" in label or " - " in label:
        strain = label.split("—")[-1].strip() if "—" in label else label.split("-")[-1].strip()
    else:
        strain = label
    marker = h.get_marker() if hasattr(h, 'get_marker') else 'o'
    strain_to_marker[strain] = marker
    print(f"  {strain:12} → {marker}")

legend_fig.legend(handles, labels, loc='center', title="Condition (color) / Strain (shape)",
                  fontsize=10, title_fontsize=11, frameon=True, ncol=6)
plt.tight_layout()
legend_fig.savefig(os.path.join(OUTPUT_DIR, "PCA_legend_separate.png"), dpi=300, bbox_inches='tight')

# ====================== 3. INDIVIDUAL SHAPES LEGEND (with extra bottom margin) ======================
print("Creating individual shapes legend...")

unique_combos = pca_df[['Strain', 'Condition']].drop_duplicates()

# Increased height and better spacing
num_rows = (len(unique_combos) + 1) // 2
fig_height = max(7, num_rows * 0.65 + 2.0)   # extra margin at bottom

fig, ax = plt.subplots(figsize=(10, fig_height))

for i, (_, row) in enumerate(unique_combos.iterrows()):
    strain = row['Strain']
    condition = row['Condition']
    color = condition_colors.get(condition, '#7f7f7f')
    marker = strain_to_marker.get(strain, 'o')   # use the extracted mapping from earlier

    col = i % 2
    row_idx = i // 2
    x = col * 1.4
    y = -row_idx * 1.25 - 0.5   # extra offset for better spacing

    ax.scatter(x, y, s=380, color=color, marker=marker, edgecolor='black', linewidth=1.5)

ax.set_xlim(-0.5, 3.0)
ax.set_ylim(-num_rows * 1.25 - 2.0, 1.5)   # extra bottom margin
ax.axis('off')

fig.savefig(os.path.join(OUTPUT_DIR, "PCA_individual_shapes.png"), dpi=300, bbox_inches='tight')
print(f"Individual shapes legend saved → PCA_individual_shapes.png")

plt.show()

print(f"\n✅ Done!")