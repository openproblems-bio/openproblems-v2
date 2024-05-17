import numpy as np


def _randomize_features(X, partition=None):
    """
    Taken and adapted from opsca-v1:
    https://github.com/openproblems-bio/openproblems/blob/acf5c95a7306b819c4a13972783433d0a48f769b/openproblems/tasks/_batch_integration/_common/methods/baseline.py#L13
    """
    X_out = X.copy()
    if partition is None:
        partition = np.full(X.shape[0], 0)
    else:
        partition = np.asarray(partition)
    for partition_name in np.unique(partition):
        partition_idx = np.argwhere(partition == partition_name).flatten()
        X_out[partition_idx] = X[np.random.permutation(partition_idx)]
    return X_out


def _randomize_graph(adata, partition=None):
    """
    Taken and adapted from opsca-v1:
    https://github.com/openproblems-bio/openproblems/blob/acf5c95a7306b819c4a13972783433d0a48f769b/openproblems/tasks/_batch_integration/_common/methods/baseline.py#L25
    """
    distances, connectivities = (
        adata.obsp["distances"],
        adata.obsp["connectivities"],
    )
    new_idx = _randomize_features(np.arange(distances.shape[0]), partition=partition)
    adata.obsp["distances"] = distances[new_idx][:, new_idx]
    adata.obsp["connectivities"] = connectivities[new_idx][:, new_idx]
    return adata


def _random_embedding(partition, jitter=0.01):
    """
    Taken and adapted from opsca-v1:
    https://github.com/openproblems-bio/openproblems/blob/acf5c95a7306b819c4a13972783433d0a48f769b/openproblems/tasks/_batch_integration/_common/methods/baseline.py#L37
    """
    from sklearn.preprocessing import LabelEncoder
    from sklearn.preprocessing import OneHotEncoder

    embedding = OneHotEncoder().fit_transform(
        LabelEncoder().fit_transform(partition)[:, None]
    )
    if jitter is not None:
        embedding = embedding + np.random.uniform(-1 * jitter, jitter, embedding.shape)
    return embedding
