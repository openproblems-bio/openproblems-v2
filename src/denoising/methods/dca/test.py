import h5py
import numpy as np
import anndata as ad
import subprocess
from os import path

input_train_path = meta["resources_dir"] + "denoising/pancreas_split_data_output_train.h5ad"
output_path = "output.h5ad"

cmd = [
    meta['executable'],
    "--input_train", input_train_path,
    "--output", output_path
]

print(">> Running script as test")
out = subprocess.run(cmd, check=True, capture_output=True, text=True)

print(">> Checking whether output file exists")
assert path.exists(output_path)

print(">> Reading h5ad files")
with h5py.File(input_train_path, 'r') as input_train:
    with h5py.File(output_path, 'r') as output:
        print("input_train:" , input_train.keys())
        print("input_train:", input_train.get("layers/counts"))
        print("output:" , output.keys())
        print("output:", output.get("layers/counts"))


        print(">> Checking whether predictions were added")
        assert "denoised" in output["layers"].keys()
        assert meta['functionality_name'] == np.string_(output["uns/method_id"]).decode('utf-8')

        print("Checking whether data from input was copied properly to output")
        assert input_train.get("layers/counts").shape == output.get("layers/counts").shape
        assert input_train["uns/dataset_id"] == output["uns/dataset_id"]

print("All checks succeeded!")
