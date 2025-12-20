"""
Functions for adding random edges to graphs for evaluation purposes
"""
import numpy as np
import torch
import pandas as pd


def add_random_edges(adata, batch_key='batch', to_graph="integration_edge_index", k_add=10, seed=None, device='cpu'):
    """
    Add random cross-batch edges to an existing graph.
    
    This function randomly selects nodes from different batches and connects them,
    without using any cell embedding information. The random edges are added to
    the existing graph stored in adata.uns[to_graph].
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    batch_key : str
        Key in adata.obs for batch labels (default: 'batch')
    to_graph : str
        Key in adata.uns for the graph edge index to modify (default: 'integration_edge_index')
    k_add : int
        Number of random cross-batch edges to add per node (default: 10)
    device : str
        Device to run on ('cpu' or 'cuda') - not used but kept for consistency
    
    Returns:
    --------
    None (modifies adata.uns[to_graph] in place)
    """
    
    # Get existing edge index
    existing_edge_index = adata.uns[to_graph]  # numpy array, shape: [2, N_edges]
    
    # Convert to numpy if it's a torch tensor
    if isinstance(existing_edge_index, torch.Tensor):
        existing_edge_index = existing_edge_index.cpu().numpy()
    existing_edge_index = np.array(existing_edge_index)
    
    # Get batch labels
    batch_labels = adata.obs[batch_key]
    if hasattr(batch_labels, 'values'):
        batch_labels = batch_labels.values
    
    # Create one-hot encoded batch matrix
    node_batch_mt = pd.get_dummies(adata.obs[batch_key]).to_numpy()
    
    N_nodes = len(batch_labels)
    
    # Get unique batches
    unique_batches = np.unique(batch_labels)
    N_batches = len(unique_batches)
    
    # Build a mapping from batch label to node indices
    batch_to_nodes = {batch: np.where(batch_labels == batch)[0] for batch in unique_batches}
    
    # Create set of existing edges for fast lookup (to avoid duplicates)
    if existing_edge_index.shape[1] > 0:
        existing_edges_set = set(map(tuple, existing_edge_index.T))
    else:
        existing_edges_set = set()
    
    # Generate random cross-batch edges
    random_edges_list = []

    if seed is not None:
         np.random.seed(seed)  # For reproducibility, can be removed if desired
    
    for node_i in range(N_nodes):
        batch_i = batch_labels[node_i]
        
        # Get all nodes from different batches
        other_batch_nodes = []
        for batch_j in unique_batches:
            if batch_j != batch_i:
                other_batch_nodes.extend(batch_to_nodes[batch_j])
        
        other_batch_nodes = np.array(other_batch_nodes)
        
        if len(other_batch_nodes) == 0:
            continue  # Skip if no other batches exist
        
        # Randomly sample k_add nodes from other batches
        n_sample = min(k_add, len(other_batch_nodes))
        sampled_targets = np.random.choice(other_batch_nodes, size=n_sample, replace=False)
        
        # Create edges and filter out duplicates
        for target_j in sampled_targets:
            edge = (node_i, target_j)
            # Only add if edge doesn't already exist (including reverse direction)
            if edge not in existing_edges_set and (target_j, node_i) not in existing_edges_set:
                random_edges_list.append(edge)
                existing_edges_set.add(edge)  # Add to set to prevent duplicates in current batch
    
    # Convert to edge index format
    if len(random_edges_list) > 0:
        random_edges = np.array(random_edges_list).T  # shape: [2, N_new_edges]
        
        # Concatenate with existing edges
        if existing_edge_index.shape[1] > 0:
            updated_edge_index = np.concatenate([existing_edge_index, random_edges], axis=1)
        else:
            updated_edge_index = random_edges
        
        # Update the graph in adata
        adata.uns[to_graph] = updated_edge_index
    else:
        # No new edges added, keep existing graph
        pass

