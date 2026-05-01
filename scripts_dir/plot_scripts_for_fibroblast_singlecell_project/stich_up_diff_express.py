import os
import pandas as pd

base_dir = os.getcwd()
frames = []

for diagnosis in os.listdir(base_dir):
    diag_path = os.path.join(base_dir, diagnosis)
    if not os.path.isdir(diag_path):
        continue
    for fname in os.listdir(diag_path):
        if not fname.endswith(".csv"):
            continue
        if fname.startswith("cell_counts_") or fname.startswith("aggregate_"):
            continue
        # e.g. "Adipogenic_vs_all.csv" -> "Adipogenic"
        cell_type = fname.replace("_vs_all.csv", "")
        df = pd.read_csv(os.path.join(diag_path, fname), index_col=0)
        df.index.name = "gene"
        df = df.reset_index()
        df.insert(0, "diagnosis", diagnosis)
        df.insert(1, "cell_type", cell_type)
        frames.append(df)

all_df = pd.concat(frames, ignore_index=True)
out_path = os.path.join(base_dir, "all.csv")
all_df.to_csv(out_path, index=False)
print(f"Written {len(all_df)} rows to {out_path}")
