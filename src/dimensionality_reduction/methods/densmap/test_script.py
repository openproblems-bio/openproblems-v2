import subprocess
import scanpy as sc
from os import path

INPUT = "toy_log_cpm_hvg.h5ad"
OUTPUT = "output.densmap.h5ad"

print(">> Running script as test")
out = subprocess.check_output([
    "./" + meta["functionality_name"],
    "--input", INPUT,
    "--output", OUTPUT
]).decode("utf-8")

print(">> Checking if output file exists")
assert path.exists(OUTPUT)

print(">> Checking if normalization method occurs")
adata = sc.read_h5ad(OUTPUT)
assert "densmap" == adata.uns["method_id"]
assert "method_code_version" in adata.uns

print(">> Running script as test")
out = subprocess.check_output([
    "./" + meta["functionality_name"],
    "--pca",
    "--input", INPUT,
    "--output", OUTPUT
]).decode("utf-8")

print(">> Checking if output file exists")
assert path.exists(OUTPUT)

print(">> Checking if normalization method occurs")
adata = sc.read_h5ad(OUTPUT)
assert "densmap" == adata.uns["method_id"]
assert "method_code_version" in adata.uns
