import anndata as ad
import sklearn.neighbors
import sklearn.pipeline
import sklearn.preprocessing
import sklearn.decomposition
import scipy.sparse

## VIASH START
par = {
    'input_train': 'resources/label_projection/openproblems_v1/pancreas.split_dataset.output_train.h5ad',
    'input_test': 'resources/label_projection/openproblems_v1/pancreas.split_dataset.output_test.h5ad',
    'output': 'output.h5ad',
    'layer_input': 'counts'
}
meta = {
    'functionality_name': 'foo',
}
## VIASH END

print("Load input data")
input_train = ad.read_h5ad(par['input_train'])
input_test = ad.read_h5ad(par['input_test'])
input_layer = par["layer_input"]

print("Set up classifier pipeline")
def pca_op(adata_train, adata_test, n_components=100):
    is_sparse = scipy.sparse.issparse(adata_train.layers[input_layer])

    min_components = min(
        [adata_train.n_obs, adata_test.n_obs, adata_train.n_vars]
    )
    if is_sparse:
        min_components -= 1
    n_components = min([n_components, min_components])
    if is_sparse:
        pca_op = sklearn.decomposition.TruncatedSVD
    else:
        pca_op = sklearn.decomposition.PCA
    return pca_op(n_components=n_components)

classifier = sklearn.neighbors.KNeighborsClassifier()
pipeline = sklearn.pipeline.Pipeline(
    [
        ("pca", pca_op(input_train, input_test, n_components=100)),
        ("scaler", sklearn.preprocessing.StandardScaler(with_mean=True)),
        ("regression", classifier),
    ]
)

print("Fit to train data")
pipeline.fit(input_train.layers[input_layer], input_train.obs["label"].astype(str))

print("Predict on test data")
input_test.obs["label_pred"] = pipeline.predict(input_test.layers[input_layer])

print("Write output to file")
input_test.uns["method_id"] = meta["functionality_name"]
input_test.write_h5ad(par['output'], compression="gzip")
