import matplotlib

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import os

OUTPUT_DIR = "output/pca_output_percent"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ====================== COLORS ======================
condition_colors = {
    'EXP': '#ff7f0e',
    'ETH': '#d62728',
    '72h': '#1f77b4',
    '96h': '#0d3b66'
}

# ====================== YOUR LIST ======================
combinations = [
    ("CENPK", "EXP"), ("CENPK", "ETH"),
    ("NER", "72h"), ("NER", "96h"),
    ("NERg", "72h"), ("NERg", "96h"),
    ("NERgE", "72h"), ("NERgE", "96h"),
    ("VAL", "72h"), ("VAL", "96h"),
    ("VALg", "72h"), ("VALg", "96h"),
    ("HUMg", "EXP"),
    ("HUMgE", "EXP"), ("HUMgE", "ETH"), ("HUMgE", "72h"), ("HUMgE", "96h"),
    ("HUMgEh", "EXP"),
    ("HUMgEhE", "EXP"), ("HUMgEhE", "ETH"), ("HUMgEhE", "72h"), ("HUMgEhE", "96h"),
    ("LIMg", "72h"), ("LIMg", "96h")
]

# ====================== EXACT MARKER MAPPING FROM YOUR GROUPED LEGEND ======================
marker_dict = {
    "CENPK": 'o',
    "NER": 'x',
    "NERg": 's',
    "NERgE": '+',
    "VAL": 'D',
    "VALg": 'd',
    "HUMg": '^',
    "HUMgE": 'x',
    "HUMgEh": 'v',
    "HUMgEhE": '*',
    "LIMg": 'h'
}

# ====================== DRAW ======================
fig, ax = plt.subplots(figsize=(10, 11))

for i, (strain, condition) in enumerate(combinations):
    color = condition_colors.get(condition, '#7f7f7f')
    marker = marker_dict.get(strain, 'o')

    col = i % 2
    row = i // 2
    x = col * 1.3
    y = -row * 1.15

    ax.scatter(x, y, s=420, color=color, marker=marker,
               edgecolor='black', linewidth=1.5)

ax.set_xlim(-0.5, 2.8)
ax.set_ylim(-len(combinations) // 2 - 2, 1)
ax.axis('off')

fig.savefig(os.path.join(OUTPUT_DIR, "PCA_individual_shapes.png"), dpi=300, bbox_inches='tight')

print("✅ Saved PCA_individual_shapes.png (now consistent with grouped legend)")
plt.show()