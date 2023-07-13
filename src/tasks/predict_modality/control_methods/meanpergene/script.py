import anndata as ad
from scipy.sparse import csc_matrix
import numpy as np

# VIASH START
par = {
    "input_train_mod1": "../../../../resources_test/predict_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.train_mod1.h5ad",
    "input_test_mod1": "../../../../resources_test/predict_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.test_mod1.h5ad",
    "input_train_mod2": "../../../../resources_test/predict_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.train_mod2.h5ad",
    "output": "../../../../resources_test/predict_modality/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.prediction.h5ad",
}

meta = {
    "functionality_name": "foo"
}
# VIASH END

input_test_mod1 = ad.read_h5ad(par["input_test_mod1"])
input_train_mod2 = ad.read_h5ad(par["input_train_mod2"])


# Find the correct shape
mean = np.array(input_train_mod2.X.mean(axis=0)).flatten()
prediction = csc_matrix(np.tile(mean, (input_test_mod1.shape[0], 1)))

# Write out prediction
out = ad.AnnData(
    X=prediction,
    obs=input_test_mod1.obs,
    var=input_train_mod2.var,
    uns={
        "dataset_id": input_test_mod1.uns["dataset_id"],
        "method_id": meta["functionality_name"],
    }
)
out.write_h5ad(par["output"], compression="gzip")
