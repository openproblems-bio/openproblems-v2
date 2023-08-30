import logging
import anndata as ad
import scipy.spatial
import scipy.sparse
import numpy as np

from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import normalize

## VIASH START
# Anything within this block will be removed by `viash` and will be
# replaced with the parameters as specified in your config.vsh.yaml.
par = {
    "input_train_mod1": "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.train_mod1.h5ad",
    "input_train_mod2": "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.train_mod2.h5ad",
    "input_train_sol": "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.train_sol.h5ad",
    "input_test_mod1": "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.test_mod1.h5ad",
    "input_test_mod2": "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.test_mod2.h5ad",
    "output": "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.prediction.h5ad",
    "n_svd": 100,
    "n_neighbors" : 10,
}

meta = {
    "funtionality_name": "foo"
}
## VIASH END

logging.basicConfig(level=logging.INFO)

logging.info("Load datasets")
input_train_mod1 = ad.read_h5ad(par["input_train_mod1"])
input_train_mod2 = ad.read_h5ad(par["input_train_mod2"])
input_train_sol = ad.read_h5ad(par["input_train_sol"])
input_test_mod1 = ad.read_h5ad(par["input_test_mod1"])
input_test_mod2 = ad.read_h5ad(par["input_test_mod2"])

# This method runs PCA on each modality individually, then runs linear regression to predict mod2
# from mod1 and finally performs kNN to match modalities

# unscramble training cells
ord = np.argsort(input_train_sol.uns['pairing_ix'])
input_train_mod2 = input_train_mod2[ord, :]

# concatenate train and test data
mod1 = ad.concat(
    {
        "train": input_train_mod1,
        "test": input_test_mod1
    },
    index_unique="-",
    label="group"
)
mod2 = ad.concat(
    {
        "train": input_train_mod2,
        "test": input_test_mod2
    },
    index_unique="-",
    label="group"
)

# Create helper views to access the test data later
mod1tr = mod1[mod1.obs["group"] == "train", :]
mod2tr = mod2[mod2.obs["group"] == "train", :]

mod1te = mod1[mod1.obs["group"] == "test", :]
mod2te = mod2[mod2.obs["group"] == "test", :]

logging.info("Running PCA")
n_svd = min(par["n_svd"], mod1.n_obs, mod2.n_obs, mod1.n_vars, mod1.n_vars)

# Use TruncatedSVD for fast decomposition of the data
mod1.obsm["X_pca"] = TruncatedSVD(n_svd).fit_transform(mod1.X)
mod2.obsm["X_pca"] = TruncatedSVD(n_svd).fit_transform(mod2.X)

reg = LinearRegression()

reg.fit(mod1tr.obsm["X_pca"], mod2tr.obsm["X_pca"])
mod2te_pred = reg.predict(mod1te.obsm["X_pca"])

neighbors = NearestNeighbors(n_neighbors=np.min((mod2te.shape[0], par["n_neighbors"])), n_jobs=-1)
neighbors = neighbors.fit(mod2te_pred)

distances, indices = neighbors.kneighbors(mod2te.obsm["X_pca"])

prediction = np.zeros((mod2te.shape[0], mod2te.shape[0]))
for i, neighbors in enumerate(indices):
    prediction[i, neighbors] = 1 / distances[i]

prediction = normalize(prediction, norm="l1")

prediction = scipy.sparse.csr_matrix(prediction)

print("Write prediction output")
prediction = ad.AnnData(
    X=prediction,
    uns={
        "dataset_id": input_train_mod1.uns["dataset_id"],
        "method_id": meta["functionality_name"]
    }
)
prediction.write_h5ad(par["output"])
