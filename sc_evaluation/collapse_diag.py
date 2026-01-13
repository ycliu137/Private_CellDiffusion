import numpy as np
from sklearn.neighbors import NearestNeighbors

def intrinsic_dimension_knn(adata, use_rep='X_dif', k=10):
    """
    Levina & Bickel (2005) intrinsic dimension estimator.
    
    Parameters
    ----------
    adata : AnnData
    use_rep : str, default 'X_dif'
        The key of the embedding to use
    k : int, number of neighbors
    
    Returns
    -------
    float
        Estimated intrinsic dimension
    """
    X = adata.obsm[use_rep]
    nbrs = NearestNeighbors(n_neighbors=k + 1).fit(X)
    distances, _ = nbrs.kneighbors(X)
    distances = distances[:, 1:]  # remove self-distance

    logs = np.log(distances[:, -1][:, None] / distances[:, :-1])
    id_est = (k - 1) / np.sum(logs, axis=1)

    return float(np.mean(id_est))

from sklearn.decomposition import PCA

def intrinsic_dimension_at_variance_percentage(adata, use_rep='X_dif', variance_percentage=0.95):
    X = adata.obsm[use_rep]
    pca = PCA().fit(X)
    cumvar = np.cumsum(pca.explained_variance_ratio_)
    n_pcs = np.searchsorted(cumvar, variance_percentage) + 1
    return n_pcs


import scipy.sparse as sp
from numpy.linalg import pinv

def variance_explained_by_embedding(adata, use_rep='X_dif'):
    """
    Fraction of variance in X explained by embedding Z via linear projection.
    
    Parameters
    ----------
    adata : AnnData
    use_rep : str, default 'X_dif'
        The key of the embedding to use
    
    Returns
    -------
    float
        R^2 variance explained
    """
    X = adata.X
    Z = adata.obsm[use_rep]

    if sp.issparse(X):
        X = X.toarray()

    # Project X onto span(Z)
    ZTZ_inv = pinv(Z.T @ Z)
    X_hat = Z @ ZTZ_inv @ Z.T @ X

    residual = X - X_hat
    r2 = 1.0 - (np.linalg.norm(residual, 'fro') ** 2 /
                np.linalg.norm(X, 'fro') ** 2)
    return float(r2)

