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

input_dir = os.path.join("data", "evolution1")
output_dir = os.path.join("output", "evolution1")

os.makedirs(output_dir, exist_ok=True)

print(f"Input directory : {input_dir}")
print(f"Output directory: {output_dir}\n")

processed_count = 0

for root, dirs, files in os.walk(input_dir):
    for filename in files:
        if not filename.lower().endswith('.fcs'):
            continue

        input_path = os.path.join(root, filename)

        try:
            meta, data = fcsparser.parse(input_path, reformat_meta=True)
            rel_path = os.path.relpath(input_path, input_dir)
            print(f"Processing: {rel_path}")

            # Extract FITC-A voltage from _channels_
            TMD = 350.0
            if '_channels_' in meta:
                channels_df = meta['_channels_']
                fitc_row = channels_df[channels_df.iloc[:, 0].str.contains('FITC-A', na=False, case=False)]
                if not fitc_row.empty:
                    try:
                        TMD = float(fitc_row.iloc[0, 4])
                        print(f"      FITC-A Voltage (TMD) = {TMD}")
                    except:
                        print(f"      Could not convert voltage, using 350")
                else:
                    print(f"      FITC-A not found, using default 350")
            else:
                print(f"      _channels_ not found, using default 350")

            # Process FITC fluorescence
            factor = ((552 / 100) ** 7.28) / ((TMD / 100) ** 7.28)
            data['normalized_FITC_A'] = data['FITC-A'] * factor
            data['FITC_fluorescence'] = (data['normalized_FITC_A'] / data['FSC-A']) / FITC_background_reference - 1

            fitc_data = data['FITC_fluorescence'].dropna().values
            if len(fitc_data) <= 1:
                print(f"   Skipped (not enough data)")
                continue

            # Kernel density
            std_fitc = np.std(fitc_data)
            bw_fitc = 10 / std_fitc if std_fitc > 0 else 1.0
            bw_fitc = max(bw_fitc, 0.12)

            kde_fitc = gaussian_kde(fitc_data, bw_method=bw_fitc)

            x = np.linspace(-500, 4000, 900)  # Max x changed to 4000
            raw_kde = kde_fitc(x)
            max_raw = np.max(raw_kde)
            scale_factor = 9.0 / max_raw if max_raw > 0 else 1.0
            y_fitc = raw_kde * scale_factor

            # ====================== CREATE IMAGE (3.3 cm width × 0.8 cm height) ======================
            fig = plt.figure(figsize=(3.3 / 2.54, 0.8 / 2.54), dpi=600)
            ax = fig.add_subplot(111)

            ax.fill_between(x, 0, y_fitc, color='green', alpha=0.40)
            ax.plot(x, y_fitc, color='green', linewidth=0.8)

            ax.set_xlim(-500, 4000)  # Max x = 4000
            ax.set_ylim(0, 10)

            # Updated x-axis ticks
            ax.set_xticks([0, 1000, 2000, 3000, 4000])
            ax.set_xticklabels([])  # No tick labels

            ax.set_yticks([])
            ax.set_ylabel('')
            ax.set_xlabel('')

            # Crisp lines
            for spine in ax.spines.values():
                spine.set_linewidth(0.6)

            # Mirror folder structure
            relative_path = os.path.relpath(root, input_dir)
            output_subdir = os.path.join(output_dir, relative_path)
            os.makedirs(output_subdir, exist_ok=True)

            output_name = os.path.splitext(filename)[0] + '.png'
            output_path = os.path.join(output_subdir, output_name)

            plt.savefig(output_path, dpi=600, bbox_inches='tight', pad_inches=0.02)
            plt.close(fig)

            processed_count += 1
            print(f"   Saved → {output_name}\n")

        except Exception as e:
            print(f"   Error processing {filename}: {e}\n")

print(f"\n✅ Batch processing complete! Generated {processed_count} individual images.")