import numpy as np
import sklearn.pipeline
import sklearn.preprocessing
import scipy.sparse
import sklearn.decomposition


def pca_op(adata_train, adata_test, n_components=100):

    is_sparse = scipy.sparse.issparse(adata_train.X)

    min_components = min(
        [adata_train.shape[0], adata_test.shape[0], adata_train.shape[1]]
    )
    if is_sparse:
        min_components -= 1
    n_components = min([n_components, min_components])
    if is_sparse:
        pca_op = sklearn.decomposition.TruncatedSVD
    else:
        pca_op = sklearn.decomposition.PCA
    return pca_op(n_components=n_components)


def classifier(adata, estimator, n_pca=100, **kwargs):
    """Run a generic scikit-learn classifier."""
    adata_train = adata[adata.obs["is_train"]]
    adata_test = adata[~adata.obs["is_train"]].copy()

    classifier = sklearn.pipeline.Pipeline(
        [
            ("pca", pca_op(adata_train, adata_test, n_components=n_pca)),
            ("scaler", sklearn.preprocessing.StandardScaler(with_mean=True)),
            ("regression", estimator(**kwargs)),
        ]
    )

    # Fit to train data
    classifier.fit(adata_train.X, adata_train.obs["celltype"].astype(str))

    # Predict on test data
    adata_test.obs["celltype_pred"] = classifier.predict(adata_test.X)

    adata.obs["celltype_pred"] = [
        adata_test.obs["celltype_pred"][idx] if idx in adata_test.obs_names else np.nan
        for idx in adata.obs_names
    ]
    return adata
