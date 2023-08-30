import anndata as ad
import numpy as np
import scipy.sparse
from sklearn.preprocessing import normalize

# VIASH START
par = {
    "input_test_mod1": "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.test_mod1.h5ad",
    "input_test_mod2": "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.test_mod2.h5ad",
    "output": "resources_test/match_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.prediction.h5ad",
}

meta = {
    "functionality_name": "foo"
}
# VIASH END

print("Load datasets")
input_test_mod1 = ad.read_h5ad(par["input_test_mod1"])
input_test_mod2 = ad.read_h5ad(par["input_test_mod2"])

# determine number of values in array
num_values = min(1000, input_test_mod1.n_obs) * input_test_mod1.n_obs
indices = np.random.randint(input_test_mod1.n_obs**2, size=num_values)

mat_x = np.random.rand(num_values)
mat_i = indices % input_test_mod1.n_obs
mat_j = (indices / input_test_mod1.n_obs).astype(int)
pairing_matrix = scipy.sparse.csr_matrix(
    (mat_x, (mat_i, mat_j)),
    shape=(input_test_mod1.n_obs, input_test_mod2.n_obs)
)

# row normalise
prob_matrix = normalize(pairing_matrix, norm="l1")

# Write out prediction
prediction = ad.AnnData(
    X=prob_matrix,
    uns={
        "method_id": meta["functionality_name"],
        "dataset_id": input_test_mod1.uns["dataset_id"]
    }
)
prediction.write_h5ad(par["output"])
