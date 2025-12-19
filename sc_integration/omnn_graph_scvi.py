import torch
import numpy
import pandas as pd
import numpy as np
import scvi

from utils.utility_fn import *


def build_omnn_scvi_graph(adata, 
                            batch_key='batch', 
                            use_rep='X_fae', 
                            n_edges_per_node=50, 
                            k_mnn=50, 
                            prune=False, 
                            device='cpu'):
    
    adata.obsm['node_batch_mt'] = pd.get_dummies(adata.obs[batch_key]).to_numpy()
    
    scvi.model.SCVI.setup_anndata(adata, batch_key= batch_key)
    model = scvi.model.SCVI(adata, n_layers=2, n_latent=50, gene_likelihood="nb")
    model.train()
    adata.obsm['X_scvi'] = model.get_latent_representation()
    
    
    N_nodes = adata.obsm[use_rep].shape[0]
    
    
    if k_mnn != 0:
        
        node_batch_mt = adata.obsm['node_batch_mt']
    
        # numpy array
        edge_index_1 = get_inter_batch_knn_edge_index(node_batch_mt, adata.obsm[use_rep], k=k_mnn)

        # torch tensor
        edge_index_1 = batch_aware_limit_outgoing_edges(torch.tensor(edge_index_1, device=device), 
                                                        torch.tensor(node_batch_mt, device=device), 
                                                        max_edges=n_edges_per_node)


        #
        ## Fuzzy clustering alignment & mnn
        #
        
        # numpy array
        edge_index_2 = get_inter_batch_knn_edge_index(node_batch_mt, adata.obsm['X_scvi'], k=k_mnn)

        
        # torch tensor
        edge_index_2 = remove_common_edges(edge_index=torch.tensor(edge_index_2, device=device), 
                                           edge_index_reference=edge_index_1)

        node_outgoing_counts_batch_mt = trg_batch_aware_node_edges_count(edge_index_1, 
                                                torch.tensor(node_batch_mt, device=device))
        
        node_batch_max_edges_mt = n_edges_per_node - node_outgoing_counts_batch_mt

        # torch tensor
        edge_index_2 = batch_and_node_aware_limit_outgoing_edges(edge_index_2, 
                                                                 torch.tensor(node_batch_mt, device=device), 
                                                                 node_batch_max_edges_mt)



        # torch tensor
        edge_index_12 = torch.cat((edge_index_1, edge_index_2), dim=1)
        edge_index_12 = limit_outgoing_edges(edge_index_12, max_edges = n_edges_per_node)
        
        
        #
        ## knn
        #

        # numpy array
        edge_index_3 = build_knn_edge_index_ckdtree(adata.obsm['X_scvi'], k=n_edges_per_node)
        
        # torch tensor
        edge_index_3 = remove_common_edges(edge_index=torch.tensor(edge_index_3, device=device), 
                                           edge_index_reference=edge_index_12)


        outgoing_counts_12, incoming_counts_12 = node_edges_count(edge_index_12, N_nodes)
        node_aware_max_edges_12 = n_edges_per_node - outgoing_counts_12
        
        # torch tensor
        edge_index_3 = node_aware_limit_outgoing_edges(edge_index_3, 
                                                       node_aware_max_edges_12)
        
        # torch tensor
        edge_index = torch.cat((edge_index_12, edge_index_3), dim=1)
    
    else:
        # numpy array
        edge_index = build_knn_edge_index_ckdtree(adata.obsm['X_scvi'], k=n_edges_per_node)
        # torch tensor
        edge_index = torch.tensor(edge_index, device=device)
        
    
    if prune:
        edge_index = prune_edges_with_IF_labels(edge_index, adata.obs['isolation'])

    
    outgoing_counts, incoming_counts = node_edges_count(edge_index, N_nodes)
    
    
    adata.uns['integration_edge_index'] = edge_index.cpu().numpy()
    
    adata.obs['incoming_counts'] = incoming_counts.cpu().numpy()
    adata.obs['outgoing_counts'] = outgoing_counts.cpu().numpy()          


def build_integration_loss_adj(adata, use_rep='X_fae', k=50, device='cpu'):
    
    node_batch_mt = extract_data_matrix_from_adata(adata, 
                                                    use_rep='node_batch_mt', 
                                                    torch_tensor=True, 
                                                    device=device)
    
    feature_matrix = extract_data_matrix_from_adata(adata, 
                                                    use_rep=use_rep, 
                                                    torch_tensor=True, 
                                                    device=device)
    N_nodes = feature_matrix.shape[0]
    
    this_array = []
    self_batch_edge_index_dict = {}
    for ii in range(node_batch_mt.shape[1]):
        batch_label = node_batch_mt[:, ii]
        batch_nodes = torch.where(batch_label == 1)[0]
        
        batch_features = feature_matrix[batch_nodes, :]
        knn_indices = feature_to_knn_indices(batch_features, 
                                             k_min=None, 
                                             k_max=k,
                                             self_included=True)
        
        edge_index = knn_indices_to_edge_index(knn_indices)
        
        label = 'batch_' + str(ii)
        
        self_batch_edge_index_dict[label] = edge_index.cpu().numpy()
        this_array.append(label)
        
    adata.uns['integration_loss_edge_index_dict'] = self_batch_edge_index_dict
    adata.uns['integration_loss_dict_index'] = this_array


def get_inter_batch_knn_edge_index(node_batch_mt, feature_matrix, k=50):
    
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
