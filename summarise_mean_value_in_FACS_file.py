import os
import pandas as pd
import fcsparser

# ====================== SETTINGS ======================
input_dir = os.path.join("data", "HUMgEh_evolution")
output_dir = os.path.join("output", "HUMgEh_evolution")

summary_output = os.path.join(output_dir, "channel_averages_summary.csv")

# List to store all results
summary_data = []

print("Starting channel average summary...")

for root, dirs, files in os.walk(input_dir):
    for filename in files:
        if not filename.lower().endswith('.fcs'):
            continue

        input_path = os.path.join(root, filename)
        rel_path = os.path.relpath(input_path, input_dir)

        try:
            meta, data = fcsparser.parse(input_path, reformat_meta=True)

            # Calculate mean for each channel
            channel_means = data.mean(numeric_only=True)

            # Prepare row for summary
            row = {
                'file_path': rel_path,
                'subfolder': os.path.relpath(root, input_dir),
                'filename': filename
            }

            # Add mean value for every channel
            for col in channel_means.index:
                row[f"mean_{col}"] = channel_means[col]

            summary_data.append(row)

            print(f"Processed averages: {rel_path}")

        except Exception as e:
            print(f"Error processing {rel_path}: {e}")

# Convert to DataFrame and save
if summary_data:
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(summary_output, index=False)
    print(f"\n✅ Channel averages summary saved to: {summary_output}")
    print(f"Total files processed: {len(summary_data)}")
else:
    print("No FCS files were processed.")