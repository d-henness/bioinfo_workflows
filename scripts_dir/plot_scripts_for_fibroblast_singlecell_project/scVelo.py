import scanpy as sc
import scvelo as scv
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filenm')
args = parser.parse_args()


adata = sc.read(args.filenm, cache=True)
print(adata.layers['spliced'])
#print(adata.layers['unspliced'])
