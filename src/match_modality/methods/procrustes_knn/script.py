import logging
import anndata as ad
import scipy.spatial
import scipy.sparse
import numpy as np

from sklearn.preprocessing import normalize
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors

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
}

meta = {
    "functionality_name": "foo"
}
## VIASH END

logging.basicConfig(level=logging.INFO)

logging.info("Load datasets")
input_train_mod1 = ad.read_h5ad(par["input_train_mod1"])
input_train_mod2 = ad.read_h5ad(par["input_train_mod2"])
# input_train_sol = ad.read_h5ad(par["input_train_sol"])
input_test_mod1 = ad.read_h5ad(par["input_test_mod1"])
input_test_mod2 = ad.read_h5ad(par["input_test_mod2"])

# This method runs PCA on each modality individually, then uses the Procrustes method to identify
# a linear transform that best superimposes the points from modality 1 onto modality 2.

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
mod1te = mod1[mod1.obs["group"] == "test", :]
mod2te = mod2[mod2.obs["group"] == "test", :]

logging.info("Running PCA")
n_svd = min(par["n_svd"], mod1.n_obs, mod2.n_obs, mod1.n_vars, mod1.n_vars)

# Use TruncatedSVD for fast decomposition of the data
mod1.obsm["X_pca"] = TruncatedSVD(n_svd).fit_transform(mod1.X)
mod2.obsm["X_pca"] = TruncatedSVD(n_svd).fit_transform(mod2.X)

logging.info("Running Procrustes Alignment")
# This function takes in two matrices of points A and B, standardizes both, and applies a linear to
# matrix B to minimize the disparity measured as the sum of the squares of the pointwise distances
# between the two input datasets
mod1.obsm["X_pro"], mod2.obsm["X_pro"], disparity = scipy.spatial.procrustes(
    mod1.obsm["X_pca"],
    mod2.obsm["X_pca"],
)
logging.info("> Disparity value is: %0.3f" % disparity)

logging.info("Perform nearest neighbors")
# To get the matching matrix, for each point in mod1_test, we take the 1000 nearest neighbors of that
# point in the transformed mod2_test dataset
n_neighbors = min(1000, mod1te.n_obs, mod1te.n_vars, mod2te.n_obs, mod2te.n_vars)
nn = NearestNeighbors(n_neighbors=n_neighbors).fit(mod1te.obsm["X_pro"])
distances, indices = nn.kneighbors(X=mod2te.obsm["X_pro"])

logging.info("Create pairing matrix")
# Translate the neighborhood assignments to a pairing matrix that is (n_obs, n_obs)
# NOTE: `pairing_matrix` must have NO MORE than 1000*n_obs non-zero entries for fast metric computation
ind_i = np.tile(np.arange(mod1te.n_obs), (n_neighbors, 1)).T.flatten()
ind_j = indices.flatten()
ind_dist = distances.flatten()
ind_x = 2 * max(ind_dist) - ind_dist
pairing_matrix = scipy.sparse.csr_matrix(
    (ind_x, (ind_i, ind_j)),
    shape=(input_test_mod1.n_obs, input_test_mod2.n_obs)
)

# row normalise
prob_matrix = normalize(pairing_matrix, norm="l1")

print("Write prediction output")
prediction = ad.AnnData(
    X=prob_matrix,
    uns={
        "dataset_id": input_train_mod1.uns["dataset_id"],
        "method_id": meta["functionality_name"]
    }
)
prediction.write_h5ad(par["output"])
