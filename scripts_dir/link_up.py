import pandas as pd
import argparse
from pathlib import Path


parser = argparse.ArgumentParser()
parser.add_argument('filenms', nargs = '+')
args = parser.parse_args()


dataframe_list = []
for filenm in args.filenms:
    new_df = pd.read_csv(filenm)
    new_df['comparison'] = Path(filenm).stem
    new_df = new_df.rename(columns = {"Unnamed: 0": "Gene"})
    dataframe_list.append(new_df)

linked_up = pd.concat(dataframe_list, ignore_index=True)
linked_up.to_csv("linked_up.csv", index = False)
