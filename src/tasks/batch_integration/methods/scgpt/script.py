import numpy as np
import anndata as ad
import json
import scipy
from sklearn.model_selection import train_test_split
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.model import TransformerModel
from scgpt.utils.util import load_pretrained
import torch

## VIASH START

par = {
    "input": "resources_test/batch_integration/pancreas/dataset.h5ad",
    "model": "resources_test/batch_integration/scgpt/pretrained_model",
    "model_config": "resources_test/batch_integration/scgpt/pretrained_model/config.json",
    "model_vocab": "resources_test/batch_integration/scgpt/pretrained_model/vocab.json",
    "pad_token": "<pad>",
    "n_bins": 51,
    "output": "output.h5ad"
}

meta = {
    "functionality_name" : "scgpt",
}

## VIASH END

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


print("Load data", flush=True)
adata = ad.read_h5ad(par["input"])

print("Preprocess data", flush=True)

if par["n_hvg"]:
    print(f"Select top {par["n_hvg"]} high variable genes", flush=True)
    idx = adata.var["hvg_score"].to_numpy().argsort()[::-1][:par["n_hvg"]]
    adata = adata[:, idx].copy()

print("Cross check genes", flush=True)

pad_token = par["pad_token"]
special_tokens = [pad_token, "<cls>", "<eoc>"]

vocab = GeneVocab.from_file(par["model_vocab"])
[vocab.append_token(s) for s in special_tokens if s not in vocab]

adata.var["id_in_vocab"] = [ 1 if feature in vocab else -1 for feature in adata.var["feature_name"]]
feature_ids_in_vocab = np.array(adata.var["id_in_vocab"])
adata = adata[:, adata.var["id_in_vocab"] >= 0]


print("Binning data", flush=True)



print("Load model config", flush=True)

with open(par["model_config"], "r") as f:
    model_configs = json.load(f)

embsize = model_configs["embsize"]
nhead = model_configs["nheads"]
d_hid = model_configs["d_hid"]
nlayers = model_configs["nlayers"]
n_layers_cls = model_configs["n_layers_cls"]

print("Tokenize input data", flush=True)
all_counts = (
    adata.layers["normalized"].A
    if scipy.sparse.issparse(adata.layers["normalized"])
    else adata.layers["normalized"]
)
genes = adata.var["feature_name"].tolist()

celltypes_labels = adata.obs["label"].tolist()  # make sure count from 0
num_types = len(set(celltypes_labels))
celltypes_labels = np.array(celltypes_labels)

batch_ids = adata.obs["batch"].tolist()
num_batch_types = len(set(batch_ids))
batch_ids = np.array(batch_ids)

(
    train_data,
    valid_data,
    train_celltype_labels,
    valid_celltype_labels,
    train_batch_labels,
    valid_batch_labels,
) = train_test_split(
    all_counts, celltypes_labels, batch_ids, test_size=0.1, shuffle=True
)

vocab.set_default_index(vocab["<pad>"])
gene_ids = np.array(vocab(genes), dtype=int)