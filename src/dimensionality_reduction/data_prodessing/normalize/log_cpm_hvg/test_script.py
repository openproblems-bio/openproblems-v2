import subprocess
import scanpy as sc
from os import path

INPUT = "toy_data.h5ad"
OUTPUT = "output.logcpmhvg.h5ad"

print(">> Running script as test")
out = subprocess.check_output([
    "./" + meta["functionality_name"],
    "--input", INPUT,
    "--n_genes", "1000",
    "--output", OUTPUT
]).decode("utf-8")

print(">> Checking if output file exists")
assert path.exists(OUTPUT)

print(">> Checking if normalization method occurs")
adata = sc.read_h5ad(OUTPUT)
assert "log_cpm_hvg" == adata.uns["normalization_method"]
