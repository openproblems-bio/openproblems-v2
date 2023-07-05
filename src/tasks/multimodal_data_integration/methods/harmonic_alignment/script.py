print("Loading dependencies")
import scanpy as sc
import harmonicalignment
import sklearn.decomposition

## VIASH START
par = {
  "mod1" : "resources_test/common/multimodal/dataset_mod1.h5ad",
  "mod2" : "resources_test/common/multimodal/dataset_mod2.h5ad",
  "output" : "output.scot.h5ad",
  "n_svd" : 100,
  "n_pca_XY" : 100,
  "eigenvectors" : 100
}
## VIASH END


print("Reading input h5ad file")
adata_mod1 = sc.read_h5ad(par["mod1"])
adata_mod2 = sc.read_h5ad(par["mod2"])

print("Check parameters")
n_svd = min([par["n_svd"], min(adata_mod1.layers["normalized"].shape) - 1, min(adata_mod2.layers["normalized"].shape) - 1])
n_eigenvectors = par["n_eigenvectors"]
n_pca_XY = par["n_pca_XY"]

if adata_mod1.layers["normalized"].shape[0] <= n_eigenvectors:
    n_eigenvectors = None
if adata_mod1.layers["normalized"].shape[0] <= n_pca_XY:
    n_pca_XY = None


print("Performing PCA reduction")
mod1_pca = sklearn.decomposition.TruncatedSVD(n_svd).fit_transform(adata_mod1.layers["normalized"])
mod2_pca = sklearn.decomposition.TruncatedSVD(n_svd).fit_transform(adata_mod2.layers["normalized"])

print("Running Harmonic Alignment")
ha_op = harmonicalignment.HarmonicAlignment(
    n_filters=8, n_pca_XY=n_pca_XY, n_eigenvectors=n_eigenvectors
)
ha_op.align(mod1_pca, mod2_pca)
XY_aligned = ha_op.diffusion_map(n_eigenvectors=n_eigenvectors)

print("Storing output data structures")

adata= sc.AnnData(
    uns={
        
    }
)

adata.obsm["aligned"] = XY_aligned[: X_pca.shape[0]]
adata.obsm["mode2_aligned"] = XY_aligned[X_pca.shape[0] :]

print("Write output to file")
adata.uns["method_id"] = "harmonic_alignment"
adata.write_h5ad(par["output"], compression = "gzip")
