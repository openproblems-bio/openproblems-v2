import anndata as ad
import numpy as np
import scipy
from sklearn.neighbors import NearestNeighbors

# VIASH START
par = {
    "input_prediction": "resources_test/joint_embedding/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.prediction.h5ad",
    "input_solution": "resources_test/joint_embedding/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.solution.h5ad",
    "output": "resources_test/joint_embedding/openproblems_bmmc_multiome_starter/openproblems_bmmc_multiome_starter.scores_totalvi.h5ad",
    "n_neighbors": 100
}
# VIASH END

print("Read input files")
predict_adata = ad.read_h5ad(par["input_prediction"])
solution_adata = ad.read_h5ad(par["input_solution"])

print("Merge prediction with solution")
merged_adata = predict_adata.copy()

batch_val = solution_adata.obs["batch"].astype(str)
batch_unique_values, batch_index = np.unique(batch_val, return_inverse=True)

merged_adata.obs["batch"] = batch_index

def entropy_batch_mixing(
    latent_space, batches, n_neighbors=50, n_pools=50, n_samples_per_pool=100
):

    def neg_kl(hist_data, global_freq):
        n_batches = len(np.unique(hist_data))
        if n_batches > 2:
            raise ValueError("Should be only two clusters for this metric")
        frequency = np.mean(hist_data == 1)
        if frequency == 0 or frequency == 1:
            return 0
        return -(
            frequency * np.log(frequency / global_freq)
            + (1 - frequency) * np.log((1 - frequency) / (1 - global_freq))
        )

    n_neighbors = min(n_neighbors, latent_space.getnnz() - 1)
    nne = NearestNeighbors(n_neighbors=1 + n_neighbors, n_jobs=8)
    nne.fit(latent_space)
    kmatrix = nne.kneighbors_graph(latent_space) - scipy.sparse.identity(
        latent_space.shape[0]
    )

    global_freq = np.mean(batches)
    print(global_freq)
    score = 0
    for t in range(n_pools):
        indices = np.random.choice(
            np.arange(latent_space.shape[0]), size=n_samples_per_pool
        )
        score += np.mean(
            [
                neg_kl(
                    batches[  # the batches of cell i's neighbors
                        kmatrix[indices].nonzero()[
                            1
                        ][  # the neighbors of cell i (columns in row i)
                            kmatrix[indices].nonzero()[0] == i  # the row of cell i
                        ]
                    ],
                    global_freq,
                )
                for i in range(n_samples_per_pool)
            ]
        )
    return score / float(n_pools)


print("Calculate latent mixing metric")
latent_mixing = entropy_batch_mixing(
    latent_space=merged_adata.obsm['X_emb'],
    batches=merged_adata.obs["batch"].values,
    n_neighbors=par["n_neighbors"]
)

print("Write output")
adata_out = ad.AnnData(
    uns = {
        "dataset_id": predict_adata.uns["dataset_id"],
        "method_id" : predict_adata.uns["method_id"],
        "metric_ids" : ["latent_mixing"],
        "metric_values" : [latent_mixing]
    }
)

adata_out.write_h5ad(par['output'], compression = "gzip")