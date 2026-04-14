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

input_dir = os.path.join("data", "fedbatch1")
output_dir = os.path.join("output", "fedbatch1")

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

            # === Check FCS version ===
            fcs_version = "Unknown"
            if '__header__' in meta:
                header = meta['__header__']
                if isinstance(header, dict) and 'FCS format' in header:
                    fcs_version = header['FCS format']
                    if isinstance(fcs_version, bytes):
                        fcs_version = fcs_version.decode('utf-8', errors='ignore')
            print(f"      FCS version: {fcs_version}")

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

            # === Save raw FITC-A and processed fluorescence (before log10) ===
            base_name = os.path.splitext(filename)[0]
            output_subdir = os.path.join(output_dir, os.path.relpath(root, input_dir))
            os.makedirs(output_subdir, exist_ok=True)

            # 1. Raw FITC-A values
            pd.DataFrame({'FITC-A': data['FITC-A']}).to_csv(
                os.path.join(output_subdir, base_name + '_FITC-A_raw.csv'), index=False)

            # 2. Processed fluorescence (before log10)
            pd.DataFrame({'FITC_fluorescence': fitc_data}).to_csv(
                os.path.join(output_subdir, base_name + '_fluorescence.csv'), index=False)

            print(f"   Saved raw FITC-A → {base_name}_FITC-A_raw.csv")
            print(f"   Saved fluorescence → {base_name}_fluorescence.csv")

            # === Log10 transformation for plotting ===
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

            x_log = np.linspace(-1.0, 5.0, 1000)
            raw_kde = kde_fitc(x_log)
            max_raw = np.max(raw_kde)
            scale_factor = 9.0 / max_raw if max_raw > 0 else 1.0
            y_fitc = raw_kde * scale_factor

            # ====================== CREATE IMAGE ======================
            fig = plt.figure(figsize=(3.3 / 2.54, 0.8 / 2.54), dpi=600)
            ax = fig.add_subplot(111)

            ax.fill_between(x_log, 0, y_fitc, color='green', alpha=0.40)
            ax.plot(x_log, y_fitc, color='green', linewidth=0.8)

            ax.set_xlim(-1.0, 5.0)
            ax.set_ylim(0, 10)

            ax.set_xticks([0, 1, 2, 3, 4, 5])
            ax.set_xticklabels([])

            ax.set_yticks([])
            ax.set_ylabel('')
            ax.set_xlabel('')

            for spine in ax.spines.values():
                spine.set_linewidth(0.6)

            # Save plot
            output_name = base_name + '.png'
            output_path = os.path.join(output_subdir, output_name)
            plt.savefig(output_path, dpi=600, bbox_inches='tight', pad_inches=0.02)
            plt.close(fig)

            processed_count += 1
            print(f"   Saved plot → {output_name}\n")

        except Exception as e:
            print(f"   Error processing {filename}: {e}\n")

print(f"\n✅ Batch processing complete! Generated {processed_count} individual images.")