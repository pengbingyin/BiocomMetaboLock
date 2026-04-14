import matplotlib

matplotlib.use('TkAgg')

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from pathlib import Path

# ====================== CONFIG ======================
INPUT_FILE = Path("data") / "Figure 5.xlsx"
OUTPUT_DIR = Path("output")
OUTPUT_DIR.mkdir(exist_ok=True)

OUTPUT_FILE = OUTPUT_DIR / "Figure_5_circle_matrix_increased_vertical.png"

WIDTH_CM = 5.0
HEIGHT_CM = 5.4  # Increased height to accommodate more vertical spacing

DPI = 600

cmap = LinearSegmentedColormap.from_list('black_to_green', ['black', '#00C853'])

# ====================== LOAD DATA ======================
df = pd.read_excel(INPUT_FILE, sheet_name="Sheet1", header=None)
values = df.iloc[:, 1:].values.astype(float)  # 20 rows × 24 columns

n_rows, n_cols = values.shape

# ====================== GRID + RADIUS ======================
radius_100_cm = 0.20  # Radius = 0.2 cm when value = 100
radius_100_in = radius_100_cm / 2.54

# Increased vertical spacing
col_spacing_cm = 0.46
row_spacing_cm = 0.54  # ← Increased vertical spacing

col_spacing = col_spacing_cm / 2.54
row_spacing = row_spacing_cm / 2.54

# ====================== PLOT ======================
fig, ax = plt.subplots(figsize=(WIDTH_CM / 2.54, HEIGHT_CM / 2.54), dpi=DPI)

for i in range(n_rows):
    for j in range(n_cols):
        val = np.clip(values[i, j], 0, 500)

        radius = min(val / 100.0, 1.0) * radius_100_in
        color = cmap(min(val / 100.0, 1.0))

        circle = plt.Circle(
            xy=(j * col_spacing, (n_rows - 1 - i) * row_spacing),  # First row at top
            radius=radius,
            color=color,
            linewidth=0.18,
            edgecolor='white'
        )
        ax.add_patch(circle)

# Clean matrix - no labels
ax.set_xticks([])
ax.set_yticks([])
ax.set_xlabel('')
ax.set_ylabel('')

ax.set_xlim(-0.15, n_cols * col_spacing + 0.1)
ax.set_ylim(-0.15, n_rows * row_spacing + 0.15)

for spine in ax.spines.values():
    spine.set_visible(False)

plt.tight_layout(pad=0.05)

plt.savefig(OUTPUT_FILE, dpi=DPI, bbox_inches='tight', pad_inches=0.02)
print(f"✅ Increased vertical spacing version saved: {OUTPUT_FILE.resolve()}")

plt.show()