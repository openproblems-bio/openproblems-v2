import subprocess
import scanpy as sc
from os import path

INPUT = "toy_data.h5ad"
OUTPUT = "result.h5ad"
MAX_ITER = "100"
NOT_REQUIRED = "test"

# Test without the 'not_required' parameter
out = subprocess.check_output(
    [
        "./" + meta['functionality_name'], # the meta is provided by Viash component
        "--input", INPUT,
        "--max_iter", MAX_ITER,
        "--output", OUTPUT
    ]
).decode("utf-8")

print(">> Checking if output file exists")
assert path.exists(OUTPUT)

print(">> Checking if data is transformed correctly")
adata = sc.read(OUTPUT)
assert "iteractions" in adata.uns
assert adata.uns["iteractions"] == 100
assert "not_required" not in adata.uns

# Test with the 'not_required' parameter
out = subprocess.check_output(
    [
        "./" + meta['functionality_name'], # the meta is provided by Viash component
        "--input", INPUT,
        "--max_iter", MAX_ITER,
        "--not_required", NOT_REQUIRED,
        "--output", OUTPUT
    ]
).decode("utf-8")

print(">> Checking if output file exists")
assert path.exists(OUTPUT)

print(">> Checking if data is transformed correctly")
adata = sc.read(OUTPUT)
assert "iteractions" in adata.uns
assert adata.uns["iteractions"] == 100
assert "not_required" in adata.uns
assert adata.uns["not_required"] == "test"
