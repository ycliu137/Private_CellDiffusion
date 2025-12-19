import torch
import numpy
import pandas as pd
import numpy as np
import scvi

from utils.utility_fn import *


def build_scvi_mnn_graph(adata, 
                    batch_key='batch', 
                    use_rep='X_fae', 
                    n_edges_per_node=50, 
                    k_mnn=50, 
                    prune=False, 
                    device='cpu'):
    """
    Build MNN (Mutual Nearest Neighbors) graph for data integration
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    batch_key : str
        Key in adata.obs for batch labels
    use_rep : str
        Key in adata.obsm for node feature representation
    n_edges_per_node : int
        Maximum number of outgoing edges per node
    k_mnn : int
        Number of mutual nearest neighbors for inter-batch connections
    prune : bool
        Whether to prune edges based on isolation forest labels
    device : str
        Device to run on ('cpu' or 'cuda')
    """

    scvi.model.SCVI.setup_anndata(adata, batch_key= batch_key)
    model = scvi.model.SCVI(adata, n_layers=2, n_latent=50, gene_likelihood="nb")
    model.train()
    adata.obsm['X_scvi'] = model.get_latent_representation()
    
    
    adata.obsm['node_batch_mt'] = pd.get_dummies(adata.obs[batch_key]).to_numpy()
    
    N_nodes = adata.obsm['X_scvi'].shape[0]
    
    node_batch_mt = adata.obsm['node_batch_mt']

    # Build inter-batch MNN edges using original representation
    # numpy array
    edge_index_1 = get_inter_batch_knn_edge_index(node_batch_mt, adata.obsm['X_scvi'], k=k_mnn)

    # torch tensor
    edge_index_1 = batch_aware_limit_outgoing_edges(torch.tensor(edge_index_1, device=device), 
                                                    torch.tensor(node_batch_mt, device=device), 
                                                    max_edges=n_edges_per_node)

    # torch tensor
    edge_index_1 = limit_outgoing_edges(edge_index_1, max_edges = n_edges_per_node)
    
    # Build final KNN edges
    # numpy array
    edge_index_3 = build_knn_edge_index_ckdtree(adata.obsm['X_scvi'], k=n_edges_per_node)
    
    # torch tensor
    edge_index_3 = remove_common_edges(edge_index=torch.tensor(edge_index_3, device=device), 
                                        edge_index_reference=edge_index_1)

    outgoing_counts_1, incoming_counts_1 = node_edges_count(edge_index_1, N_nodes)
    node_aware_max_edges_1 = n_edges_per_node - outgoing_counts_1
    
    # torch tensor
    edge_index_3 = node_aware_limit_outgoing_edges(edge_index_3, 
                                                    node_aware_max_edges_1)

    # torch tensor
    edge_index = torch.cat((edge_index_1, edge_index_3), dim=1)

        
    if prune:
        edge_index = prune_edges_with_IF_labels(edge_index, adata.obs['isolation'])

    outgoing_counts, incoming_counts = node_edges_count(edge_index, N_nodes)
    
    adata.uns['integration_edge_index'] = edge_index.cpu().numpy()
    
    adata.obs['incoming_counts'] = incoming_counts.cpu().numpy()
    adata.obs['outgoing_counts'] = outgoing_counts.cpu().numpy()


def get_inter_batch_knn_edge_index(node_batch_mt, feature_matrix, k=50):
    """
    Build inter-batch mutual nearest neighbor edge index
    
    Parameters:
    -----------
    node_batch_mt : np.ndarray
        One-hot encoded batch matrix [N_nodes, N_batches]
    feature_matrix : np.ndarray
        Feature matrix [N_nodes, N_features]
    k : int
        Number of mutual nearest neighbors
    
    Returns:
    --------
    inter_batch_edge_index : np.ndarray
        Edge index [2, N_edges] for inter-batch connections
    """
    N_nodes = node_batch_mt.shape[0]
        
    inter_batch_edge_index = numpy.empty((2, 0), dtype=np.int64)

    for ii in range(node_batch_mt.shape[1]):    # Go through all batches
        batch_trg = node_batch_mt[:, ii]
        nodes_trg = numpy.where(batch_trg == 1)[0]
        features_trg = feature_matrix[nodes_trg, :]

        nodes_src = numpy.where(batch_trg == 0)[0]
        features_src = feature_matrix[nodes_src, :]

        edge_index = build_mnn_edge_index_ckdtree(features_src, features_trg, k=k)                             

        edge_index[0,:] = nodes_src[edge_index[0,:]]
        edge_index[1,:] = nodes_trg[edge_index[1,:]]

        inter_batch_edge_index = np.concatenate((inter_batch_edge_index, edge_index), axis=1)

    return inter_batch_edge_index
