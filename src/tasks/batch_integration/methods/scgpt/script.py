import numpy as np
import anndata as ad
import json
import scipy
from sklearn.model_selection import train_test_split
from scgpt.tokenizer.gene_tokenizer import GeneVocab
from scgpt.tokenizer import tokenize_and_pad_batch
from scgpt.model import TransformerModel
from scgpt.utils.util import load_pretrained
import torch

## VIASH START

par = {
    "input": "resources_test/batch_integration/pancreas/dataset.h5ad",
    "model": "resources_test/scGPT_human/best_model.pt",
    "model_config": "resources_test/scGPT_human/args.json",
    "model_vocab": "resources_test/scGPT_human/vocab.json",
    "pad_token": "<pad>",
    "max_seq_len": None,
    "pad_value": -2,
    "n_bins": 51,
    "output": "output.h5ad",
    "n_hvg": 2000,
    "batch_size": 64,
}

meta = {
    "functionality_name" : "scgpt",
}

## VIASH END

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

def _digitize(x: np.ndarray, bins: np.ndarray) -> np.ndarray:
    assert x.ndim == 1 and bins.ndim == 1
    left_digits = np.digitize(x, bins)
    right_digits = np.digitize(x, bins, right=True)
    rands = np.random.rand(len(x))  # uniform random numbers
    digits = rands * (right_digits - left_digits) + left_digits
    digits = np.ceil(digits)
    smallest_dtype = np.min_scalar_type(digits.max().astype(np.uint)) # Already checked for non-negative values
    digits = digits.astype(smallest_dtype)
    return digits


print("Load data", flush=True)
adata = ad.read_h5ad(par["input"])

print("Preprocess data", flush=True)

if par["n_hvg"]:
    print(f"Select top {par['n_hvg']} high variable genes", flush=True)
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

layer_data = adata.layers["normalized"]

binned_rows = []
bin_edges = []
for row_number in range(adata.layers["normalized"].indptr.size-1):
    row_start_index, row_end_index = layer_data.indptr[row_number], layer_data.indptr[row_number+1]
    non_zero_row = layer_data.data[row_start_index:row_end_index]
    if non_zero_row.max() == 0:
        binned_rows.append(np.zeros_like(non_zero_row, dtype=np.int8))
        bin_edges.append(np.array([0] * par["n_bins"]))
        continue
    bins = np.quantile(non_zero_row, np.linspace(0, 1, par["n_bins"] - 1))
    non_zero_digits = _digitize(non_zero_row, bins)
    assert non_zero_digits.min() >= 1
    assert non_zero_digits.max() <= par["n_bins"] - 1
    binned_rows.append(non_zero_digits)
    bin_edges.append(np.concatenate([[0], bins]))

adata.layers["binned"] = scipy.sparse.csr_matrix((np.concatenate(binned_rows, casting="same_kind"),
                                                 layer_data.indices, layer_data.indptr), shape=layer_data.shape)

adata.obsm["bin_edges"] = np.stack(bin_edges)



print("Tokenize input data", flush=True)

all_counts = (
    adata.layers["normalized"].A
    if scipy.sparse.issparse(adata.layers["normalized"])
    else adata.layers["normalized"]
)
genes = adata.var["feature_name"].tolist()

vocab.set_default_index(vocab["<pad>"])
ntokens = len(vocab)
gene_ids = np.array(vocab(genes), dtype=int)

if not par["max_seq_len"]:
    max_seq_len = adata.var.shape[0] + 1
else:
    max_seq_len = par["max_seq_len"]

tokenized_data = tokenize_and_pad_batch(
    all_counts,
    gene_ids,
    max_len=max_seq_len,
    vocab=vocab,
    pad_token=pad_token,
    pad_value=par["pad_value"],
    append_cls=True,  # append <cls> token at the beginning,
    include_zero_gene=False,
    return_pt=True,
    mod_type=None,
    vocab_mod=None
    )

all_gene_ids, all_values = tokenized_data["genes"], tokenized_data["values"]
padding_mask = all_gene_ids.eq(vocab[pad_token])
padding_mask = padding_mask.numpy()

print("Load model config", flush=True)

with open(par["model_config"], "r") as f:
    model_configs = json.load(f)

embsize = model_configs["embsize"]
nhead = model_configs["nheads"]
d_hid = model_configs["d_hid"]
nlayers = model_configs["nlayers"]
n_layers_cls = model_configs["n_layers_cls"]

batch_id_cats = adata.obs["batch"].astype("category")
batch_id_labels = batch_id_cats.cat.codes.values
batch_ids = batch_id_labels.tolist()
batch_ids = np.array(batch_ids)
num_batch_types = len(set(batch_ids))

model = TransformerModel(
    ntokens,
    d_model=embsize,
    nhead=nhead,
    d_hid=d_hid,
    nlayers=nlayers,
    vocab=vocab,
    dropout=0.5, # scGPT default, only relevant for fine-tuning applications
    pad_token=pad_token,
    pad_value=par["pad_value"],
    nlayers_cls=3,  # only applicable for decoder-based operations
    n_cls=1,  # only applicable for decoder-based operations
    do_mvc=False,  # only applicable for decoder-based operations
    ecs_threshold=0.8,  # only applicable for decoder-based operations
    do_dab=False,  # only applicable for decoder-based operations
    use_batch_labels=False, # only applicable for decoder-based operations
    num_batch_labels=num_batch_types,
    domain_spec_batchnorm=True,
    input_emb_style="continuous",  # scGPT default
    explicit_zero_prob=False,  #TODO: Parametrize when GPU-based machine types are supported
    use_fast_transformer=True if device == "cuda" else False,  #TODO: Parametrize when GPU-based machine types are supported
    # fast_transformer_backend="flash",  #TODO: Parametrize when GPU-based machine types are supported
    pre_norm=False  #TODO: Parametrize when GPU-based machine types are supported
    )

load_pretrained(
    model,
    torch.load(par["model"], map_location=device),
    verbose=False
    )

model.to(device)
model.eval()


cell_embeddings = model.encode_batch(
    torch.from_numpy(np.array(all_gene_ids)),
    torch.from_numpy(np.array(all_values)).float(),
    src_key_padding_mask=torch.from_numpy(np.array(padding_mask)),
    batch_size=par["batch_size"],
    batch_labels=torch.from_numpy(batch_ids).long(),
    output_to_cpu=True,
    time_step=0,
    return_np=True
)

cell_embeddings = cell_embeddings / np.linalg.norm(
    cell_embeddings, axis=1, keepdims=True
)

print("Store outputs", flush=True)
output = ad.AnnData(
    obs=adata.obs[[]],
    var=adata.var[[]],
    obsm={
        "X_emb": cell_embeddings,
    },
    uns={
        "dataset_id": adata.uns["dataset_id"],
        "normalization_id": adata.uns["normalization_id"],
        "method_id": meta["functionality_name"],
    },
)

print("Write output to file", flush=True)
output.write_h5ad(par["output"], compression="gzip")

# celltypes_labels = adata.obs["label"].tolist()  # make sure count from 0
# num_types = len(set(celltypes_labels))
# celltypes_labels = np.array(celltypes_labels)

# batch_ids = adata.obs["batch"].tolist()
# num_batch_types = len(set(batch_ids))
# batch_ids = np.array(batch_ids)

# (
#     train_data,
#     valid_data,
#     train_celltype_labels,
#     valid_celltype_labels,
#     train_batch_labels,
#     valid_batch_labels,
# ) = train_test_split(
#     all_counts, celltypes_labels, batch_ids, test_size=0.1, shuffle=True
# )

