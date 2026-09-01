import fcsparser
import os
import pandas as pd
import matplotlib

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import numpy as np

# ====================== SETTINGS ======================
plt.rcParams['font.family'] = 'Arial'

FITC_background_reference = 0.009643887

folder_name = "nerolidol_figure2"
input_dir = os.path.join("data", folder_name)
output_dir = os.path.join("output", folder_name)

desired_tubes = ['Tube_010', 'Tube_011', 'Tube_012']
time_points = ['2h', '6h', '10h', '24h', '48h', '72h']

os.makedirs(output_dir, exist_ok=True)

tube_data = {t: [] for t in desired_tubes}

print(f"Input directory : {input_dir}")
print(f"Output directory: {output_dir}\n")

processed_count = 0

for root, dirs, files in os.walk(input_dir):
    relative_path = os.path.relpath(root, input_dir).strip(os.sep)
    if relative_path not in time_points:
        continue

    print(f"Found time-point folder: {relative_path}")

    for filename in files:
        if not filename.lower().endswith('.fcs'):
            continue
        if not any(t in filename for t in desired_tubes):
            continue

        input_path = os.path.join(root, filename)
        try:
            meta, data = fcsparser.parse(input_path, reformat_meta=True)

            TMD = 400 if relative_path in ['2h', '6h', '10h'] else 350
            factor = ((552 / 100) ** 7.28) / ((TMD / 100) ** 7.28)
            data['normalized_FITC_A'] = data['FITC-A'] * factor

            data['FITC_fluorescence'] = (data['normalized_FITC_A'] / data['FSC-A']) / FITC_background_reference - 1

            fitc_data = data['FITC_fluorescence'].dropna().values
            if len(fitc_data) > 1:
                for t in desired_tubes:
                    if t in filename:
                        tube_data[t].append((relative_path, filename, fitc_data))
                        break

            base_name = os.path.splitext(filename)[0]
            output_subdir = os.path.join(output_dir, relative_path)
            os.makedirs(output_subdir, exist_ok=True)
            data.to_csv(os.path.join(output_subdir, base_name + '.csv'), index=False)

            processed_count += 1

        except Exception as e:
            print(f"   Error processing {filename}: {e}")

print(f"\n✅ Finished! Successfully processed {processed_count} files.\n")

# ====================== CREATE FINAL FIGURE ======================
fig_width_cm = 15.0
subplot_height_cm = 1.85
fig_height_cm = subplot_height_cm * len(time_points)

fig, axs = plt.subplots(len(time_points), 3,
                        figsize=(fig_width_cm/2.54, fig_height_cm/2.54),
                        sharex=True, sharey=True,
                        gridspec_kw={'hspace': 0.28, 'wspace': 0.15})

label_fontsize = 8
tick_fontsize = 8

for row, time_point in enumerate(time_points):
    for col, tube in enumerate(desired_tubes):
        ax = axs[row, col]

        sample = next((s for s in tube_data[tube] if s[0] == time_point), None)
        if sample is None:
            ax.axis('off')
            continue

        _, filename, fitc_data = sample

        # === Log10 transformation as requested ===
        fitc_log = np.zeros_like(fitc_data, dtype=float)
        for i, val in enumerate(fitc_data):
            if val > 0:
                fitc_log[i] = np.log10(val)
            elif val < 0:
                fitc_log[i] = -np.log10(-val)
            else:
                fitc_log[i] = 0.0

        # Kernel density on log-transformed data
        std_log = np.std(fitc_log)
        bw_log = 0.15 / std_log if std_log > 0 else 0.2
        bw_log = max(bw_log, 0.08)

        kde_fitc = gaussian_kde(fitc_log, bw_method=bw_log)

        # X-scale: min = -1, max = 5
        x_log = np.linspace(-1.0, 5.0, 1000)
        raw_kde = kde_fitc(x_log)

        max_raw = np.max(raw_kde)
        scale_factor = 9.0 / max_raw if max_raw > 0 else 1.0
        y_fitc = raw_kde * scale_factor

        ax.fill_between(x_log, 0, y_fitc, color='green', alpha=0.35)
        ax.plot(x_log, y_fitc, color='green', linewidth=1.1)

        ax.set_xlim(-1.0, 5.0)
        ax.set_ylim(0, 10)

        ax.set_yticks([])
        ax.set_ylabel('')

        if row == len(time_points) - 1:
            ax.set_xticks([0, 1, 2, 3, 4, 5])
            ax.set_xticklabels(['0', '1', '2', '3', '4', '5'],
                               fontsize=tick_fontsize,
                               rotation=30, ha='center')
            # Small right shift
            for label in ax.get_xticklabels():
                label.set_x(label.get_position()[0] + 0.015)
        else:
            ax.set_xticks([])

        # Time point label (normal, size 8)
        ax.text(0.97, 0.96, time_point,
                transform=ax.transAxes,
                fontsize=label_fontsize,
                ha='right', va='top')

plt.tight_layout()
plot_path = os.path.join(output_dir, 'combined_density_final.png')
plt.savefig(plot_path, dpi=300, bbox_inches='tight')
print(f"Saved final figure (15 cm width) → {plot_path}")

print("\nBatch processing complete.")