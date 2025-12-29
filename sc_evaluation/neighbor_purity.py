"""
Evaluation functions for single-cell data integration quality
"""
import numpy as np
import torch
import pandas as pd
from utils.utility_fn import build_mnn_edge_index_ckdtree, build_knn_edge_index_ckdtree


def evaluate_neighbor_purity(adata, label_key='labels', graph_key='integration_edge_index'):
    """
    Calculate the proportion of edges connecting cells with the same label.
    
    This metric measures how well the graph preserves cell type identity:
    higher values indicate that similar cells (same label) are connected,
    which is desirable for data integration quality.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    label_key : str
        Key in adata.obs for cell labels (default: 'labels')
    graph_key : str
        Key in adata.uns for edge index (default: 'integration_edge_index')
    
    Returns:
    --------
    purity : float
        Proportion of edges connecting cells with the same label (0.0 to 1.0)
    """
    
    edge_index = adata.uns[graph_key]  # numpy array, shape: [2, N_edges]
    labels = adata.obs[label_key]
    
    # Convert labels to numpy array if it's a pandas Series
    if hasattr(labels, 'values'):
        labels = labels.values
    
    # Convert edge_index to numpy if it's a torch tensor
    if isinstance(edge_index, torch.Tensor):
        edge_index = edge_index.cpu().numpy()
    
    # Ensure edge_index is 2D array
    edge_index = np.array(edge_index)
    if edge_index.ndim == 1:
        raise ValueError(f"edge_index should be 2D array with shape [2, N_edges], got shape {edge_index.shape}")
    
    # Get source and target node indices
    source_nodes = edge_index[0, :]  # shape: [N_edges]
    target_nodes = edge_index[1, :]  # shape: [N_edges]
    
    # Get labels for source and target nodes
    source_labels = labels[source_nodes]  # shape: [N_edges]
    target_labels = labels[target_nodes]  # shape: [N_edges]
    
    # Count edges where source and target have the same label
    same_label_mask = source_labels == target_labels
    num_same_label_edges = np.sum(same_label_mask)
    
    # Total number of edges
    total_edges = edge_index.shape[1]
    
    # Calculate purity
    if total_edges == 0:
        purity = 0.0
    else:
        purity = num_same_label_edges / total_edges
    
    return purity


def evaluate_mnn_neighbor_purity(adata, use_rep='X_fae', batch_key='batch', label_key='labels', k_mnn=50):
    """
    Build MNN graph and calculate neighbor purity based on cell labels.
    
    This function builds a mutual nearest neighbor (MNN) graph between different batches
    and then evaluates how well the graph preserves cell type identity by calculating
    the proportion of edges connecting cells with the same label.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    use_rep : str
        Key in adata.obsm for feature representation (default: 'X_fae')
    batch_key : str
        Key in adata.obs for batch labels (default: 'batch')
    label_key : str
        Key in adata.obs for cell labels (default: 'labels')
    k_mnn : int
        Number of mutual nearest neighbors for inter-batch connections (default: 50)
    
    Returns:
    --------
    purity : float
        Proportion of edges connecting cells with the same label (0.0 to 1.0)
    edge_index : np.ndarray
        The MNN edge index that was built [2, N_edges]
    """
    
    data_matrix = adata.obsm[use_rep]
    batch_labels = adata.obs[batch_key]
    labels = adata.obs[label_key]

    # Convert labels to numpy array if it's a pandas Series
    if hasattr(labels, 'values'):
        labels = labels.values
    
    # Convert batch_labels to numpy array if it's a pandas Series
    if hasattr(batch_labels, 'values'):
        batch_labels = batch_labels.values
    
    # Ensure data_matrix is numpy array
    if isinstance(data_matrix, torch.Tensor):
        data_matrix = data_matrix.cpu().numpy()
    data_matrix = np.array(data_matrix)
    
    # Create one-hot encoded batch matrix
    # Convert batch_labels to DataFrame temporarily to use pd.get_dummies
    node_batch_mt = pd.get_dummies(adata.obs[batch_key]).to_numpy()
    
    # Build inter-batch MNN edge index
    N_nodes = data_matrix.shape[0]
    inter_batch_edge_index = np.empty((2, 0), dtype=np.int64)
    
    # Go through all batches to build inter-batch MNN edges
    for ii in range(node_batch_mt.shape[1]):
        batch_trg = node_batch_mt[:, ii]
        nodes_trg = np.where(batch_trg == 1)[0]
        features_trg = data_matrix[nodes_trg, :]
        
        nodes_src = np.where(batch_trg == 0)[0]
        features_src = data_matrix[nodes_src, :]
        
        # Skip if no source or target nodes
        if len(nodes_src) == 0 or len(nodes_trg) == 0:
            continue
        
        # Build MNN edges between source and target batches
        edge_index = build_mnn_edge_index_ckdtree(features_src, features_trg, k=k_mnn)
        
        # Map back to original node indices
        edge_index[0, :] = nodes_src[edge_index[0, :]]
        edge_index[1, :] = nodes_trg[edge_index[1, :]]
        
        inter_batch_edge_index = np.concatenate((inter_batch_edge_index, edge_index), axis=1)
    
    edge_index = inter_batch_edge_index
    
    # Calculate neighbor purity using the built MNN graph
    if edge_index.shape[1] == 0:
        purity = 0.0
    else:
        # Get source and target node indices
        source_nodes = edge_index[0, :]  # shape: [N_edges]
        target_nodes = edge_index[1, :]  # shape: [N_edges]
        
        # Get labels for source and target nodes
        source_labels = labels[source_nodes]  # shape: [N_edges]
        target_labels = labels[target_nodes]  # shape: [N_edges]
        
        # Count edges where source and target have the same label
        same_label_mask = source_labels == target_labels
        num_same_label_edges = np.sum(same_label_mask)
        
        # Total number of edges
        total_edges = edge_index.shape[1]
        
        # Calculate purity
        purity = num_same_label_edges / total_edges
    
    return purity, edge_index
        

def evaluate_knn_neighbor_purity(adata, use_rep='X_fae', label_key='labels', k=50):
    """
    Calculate the proportion of edges connecting cells with the same label.
    
    This metric measures how well the graph preserves cell type identity:
    higher values indicate that similar cells (same label) are connected,
    which is desirable for data integration quality.
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    use_rep : str
        Key in adata.obsm for feature representation (default: 'X_fae')
    label_key : str
        Key in adata.obs for cell labels (default: 'labels')
    k : int
        Number of nearest neighbors for KNN graph (default: 50)
    
    Returns:
    --------
    purity : float
        Proportion of edges connecting cells with the same label (0.0 to 1.0)
    edge_index : np.ndarray
        The KNN edge index that was built [2, N_edges]
    """
    data_matrix = adata.obsm[use_rep]
    labels = adata.obs[label_key]

    # Convert labels to numpy array if it's a pandas Series
    if hasattr(labels, 'values'):
        labels = labels.values
    
    # Ensure data_matrix is numpy array
    if isinstance(data_matrix, torch.Tensor):
        data_matrix = data_matrix.cpu().numpy()
    data_matrix = np.array(data_matrix)

    # Build KNN graph
    edge_index = build_knn_edge_index_ckdtree(data_matrix, k=k)

    # Calculate neighbor purity using the built KNN graph
    if edge_index.shape[1] == 0:
        purity = 0.0
    else:
        # Get source and target node indices
        source_nodes = edge_index[0, :]  # shape: [N_edges]
        target_nodes = edge_index[1, :]  # shape: [N_edges]

        # Get labels for source and target nodes
        source_labels = labels[source_nodes]  # shape: [N_edges]
        target_labels = labels[target_nodes]  # shape: [N_edges]
        
        # Count edges where source and target have the same label
        same_label_mask = source_labels == target_labels
        num_same_label_edges = np.sum(same_label_mask)
        
        # Total number of edges
        total_edges = edge_index.shape[1]
        
        # Calculate purity
        purity = num_same_label_edges / total_edges
    
    return purity, edge_index