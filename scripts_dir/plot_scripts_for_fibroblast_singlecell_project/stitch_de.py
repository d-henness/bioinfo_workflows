#!/usr/bin/env python3
import os
import re
import pandas as pd

base_dir = os.path.dirname(os.path.abspath(__file__))
frames = []

for diagnosis in os.listdir(base_dir):
    diag_path = os.path.join(base_dir, diagnosis)
    if not os.path.isdir(diag_path):
        continue
    for fname in os.listdir(diag_path):
        if not fname.endswith('.csv'):
            continue
        # Only DE result files: {diagnosis}__vs_normal_skin_{CellType}.csv
        m = re.fullmatch(rf'{re.escape(diagnosis)}__vs_normal_skin_(.+)\.csv', fname)
        if not m:
            continue
        cell_type = m.group(1)
        df = pd.read_csv(os.path.join(diag_path, fname), index_col=0)
        df.index.name = 'gene'
        df = df.reset_index()
        df.insert(0, 'diagnosis', diagnosis)
        df.insert(1, 'cell_type', cell_type)
        frames.append(df)

result = pd.concat(frames, ignore_index=True)
out_path = os.path.join(base_dir, 'all.csv')
result.to_csv(out_path, index=False)
print(f"Wrote {len(result)} rows to {out_path}")
