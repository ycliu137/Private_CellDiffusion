"""
Graph Convolutional Network (GCN) based integration module
This is an alternative to the diffusion-based integration in integration_DIF.py
"""
import torch
import torch.nn as nn
import torch.nn.functional as F

from utils import info_log
from utils.utility_fn import extract_data_matrix_from_adata

from sc_integration.integration_graph import build_integration_graph, build_integration_loss_adj


def integration_gcn(adata, 
                    use_rep='X_fae', 
                    save_key='X_gcn', 
                    max_epoch=2000, 
                    lr=1e-3, 
                    device='cpu',
                    num_features_gcn=50,
                    num_layers_gcn=2,
                    activation=nn.ELU(),
                    data_dtype = torch.float32,
                    dropout=0.0, 
                    encoder=None, 
                    decoder=[300],
                    save_model = True,
                    load_model_state = False,
                    loss_reduction = "sum"):
    """
    GCN-based data integration function
    
    Parameters:
    -----------
    adata : AnnData
        Annotated data object
    use_rep : str
        Key for node features representation
    save_key : str
        Key to save the integrated embeddings in adata.obsm
    max_epoch : int
        Maximum number of training epochs
    lr : float
        Learning rate
    device : str
        Device to run on ('cpu' or 'cuda')
    num_features_gcn : int
        Number of hidden features in GCN layers
    num_layers_gcn : int
        Number of GCN layers
    activation : nn.Module
        Activation function
    data_dtype : torch.dtype
        Data type for tensors
    dropout : float
        Dropout rate
    encoder : list or None
        Encoder layer dimensions
    decoder : list or None
        Decoder layer dimensions
    save_model : bool
        Whether to save model state
    load_model_state : bool
        Whether to load model state from adata.uns
    loss_reduction : str
        Loss reduction method ('sum' or 'mean')
    """
    
    gcn_args = {"use_rep": use_rep,
                 "num_features_gcn": num_features_gcn,
                 "num_layers_gcn": num_layers_gcn,
                 "dropout": dropout, 
                 "encoder": encoder, 
                 "decoder": decoder,
                 "save_model": save_model,
                 "load_model_state": load_model_state}
    
    
    info_log.print('--------> Starting GCN-based data integration ...')
    
    # data
    feature_matrix = extract_data_matrix_from_adata(adata, 
                                                    use_rep=use_rep, 
                                                    torch_tensor=True, 
                                                    data_dtype=data_dtype, 
                                                    device=device)
        
    edge_index = torch.tensor(adata.uns['integration_edge_index'], dtype=torch.int64, device=device)
    
    num_of_nodes = feature_matrix.shape[0]
    
    node_batch_mt = torch.tensor(adata.obsm['node_batch_mt'], device=device)
    node_batch_mask = node_batch_mt.to(torch.bool)
    
    adjacency_list = []
    
    for ii in range(node_batch_mt.size(1)):
        N_nodes_batch = node_batch_mt[:,ii].sum().item()
        batch_dict_index = adata.uns['integration_loss_dict_index'][ii]
        edge_index_now = torch.tensor(adata.uns['integration_loss_edge_index_dict'][batch_dict_index], 
                                      dtype=torch.int64)

        adjacency_now = edge_index_to_adj(edge_index_now.to(device), N_nodes_batch)
        adjacency_list.append(adjacency_now)

    D_in = feature_matrix.shape[1]
    D_out = D_in
    
    if encoder is None:
        encoder = None if D_in==num_features_gcn else [D_in, num_features_gcn]
    else:
        encoder = [D_in] + encoder + [num_features_gcn]
    
    if decoder is None:
        decoder = None if D_out==num_features_gcn else [num_features_gcn, D_out]
    else:
        decoder = [num_features_gcn] + decoder + [D_out]
        

    model_gcn = Graph_GCN(num_features_gcn = num_features_gcn, 
                          num_layers_gcn=num_layers_gcn,
                          activation=activation,
                          dropout=dropout, 
                          encoder=encoder, 
                          decoder=decoder,
                          node_batch_one_hot=node_batch_mt).to(device)

    if load_model_state:
        try: 
            state_dict_torch = {k: torch.tensor(v).to(device) for k, v in adata.uns['gcn_state_dict'].items()}
            model_gcn.load_state_dict(state_dict_torch)
        except:
            print("GCN model failed to load model state.")
                            
    optimizer = torch.optim.Adam(model_gcn.parameters(), lr=lr)

    # Build normalized adjacency matrix for GCN
    adj_normalized = normalize_adjacency(edge_index, num_of_nodes, device=device)

    for epoch in range(max_epoch):
        model_gcn.train()
        optimizer.zero_grad()
        
        data = (feature_matrix, adj_normalized)

        out_nodes_features, recon_adj_list, last_embedding = model_gcn(data)
        
        loss_list = []
        for ii in range(len(adjacency_list)):
            target_now = torch.tensor(adjacency_list[ii].to(device), dtype = recon_adj_list[ii].dtype)
            loss_now = F.binary_cross_entropy_with_logits(recon_adj_list[ii], target_now, reduction=loss_reduction)
            loss_list.append(loss_now)
            
        loss = loss_list[0]
        for ii in range(1, len(loss_list)):
            loss += loss_list[ii]
        

        # Backprop and Update
        loss.backward()
        cur_loss = loss.item()
        optimizer.step()
        
        if epoch%50 == 0:
            info_log.interval_print(f"----------------> Epoch: {epoch+1}/{max_epoch}, Current loss: {cur_loss:.4f}")
    
    info_log.interval_print(f"----------------> Epoch: {epoch+1}/{max_epoch}, Current loss: {cur_loss:.4f}")

    # save model state
    if save_model:
        state_dict_numpy = {k: v.detach().cpu().numpy() for k, v in model_gcn.state_dict().items()}
        adata.uns['gcn_state_dict'] = state_dict_numpy

    if save_key is None:
        adata.obsm['X_gcn'] = last_embedding.detach().cpu().numpy()
    else:
        adata.obsm[save_key] = last_embedding.detach().cpu().numpy()
        
    adata.uns['gcn_args'] = gcn_args
        

class Graph_GCN(nn.Module):
    """
    Graph Convolutional Network for data integration
    """
    def __init__(self, num_features_gcn,
                 num_layers_gcn=2, 
                 activation=nn.ELU(),
                 dropout=0.0,  
                 encoder=None, 
                 decoder=None,
                 node_batch_one_hot=None):
        super().__init__()
        
        self.node_batch_mask = node_batch_one_hot.to(torch.bool)
        
        self.gcn = GCN(num_features_gcn = num_features_gcn, 
                       num_layers_gcn=num_layers_gcn,
                       activation=activation,
                       dropout=dropout, 
                       encoder=encoder, 
                       decoder=decoder)

        self.decode = InnerProductDecoder(0, act=lambda x: x)


    def forward(self, data):
        
        data, last_embedding = self.gcn(data)
        
        out_nodes_features, adj_normalized = data
        
        recon_adj_list = []
        for ii in range(self.node_batch_mask.size(1)):
            
            out_features_now = out_nodes_features[self.node_batch_mask[:, ii],:]
            recon_adj_now = self.decode(out_features_now)
            recon_adj_list.append(recon_adj_now)
        
        return out_nodes_features, recon_adj_list, last_embedding
        

class GCN(nn.Module):
    """
    Graph Convolutional Network model
    """
    def __init__(self, num_features_gcn, num_layers_gcn=2,
                 activation=nn.ELU(),
                 dropout=0.0,  
                 encoder=None, decoder=None):
        super().__init__()

        self.num_features_gcn = num_features_gcn
        self.num_layers_gcn = num_layers_gcn
        
        self.activation = activation
        self.encoder = nn.Identity()
        self.decoder = nn.Identity()
        
        # Encoder: change input dimension to GCN dimension
        if encoder is not None:
            assert encoder[-1] == self.num_features_gcn, f'"encoder" should be a vector with the length represents '\
                    f'the number of encoder layers and the components represents the dimension of features in each layer. '\
                    f'Note the first component should be the dimension of the input features and the last component '\
                    f'should equal to "num_features_gcn".'
            self.encoder = AutoEncoder(encoder, activation = self.activation, last_activation=False)
        
        # Decoder: change GCN dimension to output dimension
        if decoder is not None:
            assert decoder[0] == self.num_features_gcn, f'"decoder" should be a vector with the length represents '\
                    f'the number of decoder layers and the components represents the dimension of features in each layer. '\
                    f'Note the first component should equal to "num_features_gcn" and the last component '\
                    f'should be the dimension of the output features.'
            self.decoder = AutoEncoder(decoder, activation = self.activation, last_activation=False)

        # GCN layers
        self.gcn_layers = nn.ModuleList()
        self.dropout = nn.Dropout(p=dropout)
        
        # Create GCN layers
        for i in range(num_layers_gcn):
            self.gcn_layers.append(GCNLayer(num_features_gcn, activation=activation))

    def forward(self, data):
        """
        data is (in_nodes_features, adj_normalized) tuple
        """
        #
        # Step 1: Data pre-processing and encoding
        #
        
        # Change the node feature dimension: [N, NFIN] -> [N, NFGCN]
        data = self.encoder(data)
        
        in_nodes_features, adj_normalized = data
        
        #
        # Step 2: GCN Forward Pass
        #
        
        # Apply GCN layers
        x = in_nodes_features
        for i, gcn_layer in enumerate(self.gcn_layers):
            x = gcn_layer(x, adj_normalized)
            if i < len(self.gcn_layers) - 1:  # Apply dropout except for last layer
                x = self.dropout(x)
        
        last_embedding = x
        
        #
        # Step 3: Data decoding and outputting
        #
        
        # Change the node feature dimension: [N, NFGCN] -> [N, NFOUT]
        data = self.decoder((x, adj_normalized))
        
        return data, last_embedding


class GCNLayer(nn.Module):
    """
    Graph Convolutional Layer
    Implements: H^(l+1) = σ(D^(-1/2) A D^(-1/2) H^(l) W^(l))
    """
    def __init__(self, num_features, activation=nn.ELU()):
        super().__init__()
        self.num_features = num_features
        self.linear = nn.Linear(num_features, num_features, bias=False)
        self.activation = activation
        nn.init.xavier_uniform_(self.linear.weight)
        
    def forward(self, x, adj_normalized):
        """
        x: node features [N, F]
        adj_normalized: normalized adjacency matrix [N, N]
        """
        # Graph convolution: A_norm @ X @ W
        x = torch.matmul(adj_normalized, x)
        x = self.linear(x)
        x = self.activation(x)
        return x


class AutoEncoder(nn.Module):
    """
    Autoencoder for encoding/decoding node features
    Same as in diffusion/gnd.py
    """
    def __init__(self, num_features_list, activation, last_activation, pre_activation=False):
        super().__init__()
        self.activation = activation
        
        linear_layers = []
        if pre_activation:  # Add an activation layer in front of autoencoder
            linear_layers.append(self.activation)
            
        for i in range(len(num_features_list)-1):
            N_in = num_features_list[i]
            N_out = num_features_list[i+1]
            layer = nn.Linear(N_in, N_out, bias=False)
            
            # The default TF initialization
            nn.init.xavier_uniform_(layer.weight)
            
            linear_layers.append(layer)
            linear_layers.append(self.activation)
            
        if last_activation==False:  # Remove the last activation layer in the autoencoder
            linear_layers = linear_layers[:-1]
            
        self.autoencoder = nn.Sequential(*linear_layers,)
             
    def forward(self, data): 
        in_nodes_features, adj_normalized = data
        out_nodes_features = self.autoencoder(in_nodes_features)
        return out_nodes_features, adj_normalized


class InnerProductDecoder(nn.Module):
    """
    Inner product decoder for reconstructing adjacency matrix
    Same as in integration_DIF.py
    """
    def __init__(self, dropout, act=torch.sigmoid):
        super().__init__()
        self.dropout = dropout
        self.act = act

    def forward(self, z):
        z = F.dropout(z, self.dropout, training=self.training)
        adj = self.act(torch.mm(z, z.t()))
        return adj
    

def edge_index_to_adj(edge_index, num_of_nodes):
    """
    Construct adjacency matrix from edge index
    Same as in integration_DIF.py
    """
    adjacency_matrix = torch.zeros((num_of_nodes, num_of_nodes), dtype=edge_index.dtype, device=edge_index.device)
    adjacency_matrix[edge_index[0], edge_index[1]] = 1
    return adjacency_matrix


def normalize_adjacency(edge_index, num_of_nodes, device='cpu'):
    """
    Normalize adjacency matrix for GCN: D^(-1/2) A D^(-1/2)
    where D is the degree matrix and A is the adjacency matrix
    """
    # Create adjacency matrix
    adj = edge_index_to_adj(edge_index, num_of_nodes)
    
    # Convert to float dtype for normalization calculations
    adj = adj.float()
    
    # Add self-loops
    adj = adj + torch.eye(num_of_nodes, device=device, dtype=adj.dtype)
    
    # Calculate degree matrix
    degree = adj.sum(dim=1)
    degree_inv_sqrt = torch.pow(degree, -0.5)
    degree_inv_sqrt[degree_inv_sqrt == float('inf')] = 0  # Handle isolated nodes
    
    # Normalize: D^(-1/2) A D^(-1/2)
    degree_matrix_inv_sqrt = torch.diag(degree_inv_sqrt)
    adj_normalized = torch.matmul(torch.matmul(degree_matrix_inv_sqrt, adj), degree_matrix_inv_sqrt)
    
    return adj_normalized


def integration_gcn_high_throughput_mode(adata, 
                                         batch_key='batch', 
                                         use_rep='X_fae', 
                                         max_epoch=2000, 
                                         lr=1e-3, 
                                         device='cpu'):
    """
    High-throughput mode for GCN-based integration
    Similar to integration_high_throughput_mode in integration_DIF.py
    """
    
    info_log.print('--------> Build integration graph for GCN ...')
    
    build_integration_graph(adata, 
                            batch_key=batch_key, 
                            use_rep=use_rep, 
                            n_edges_per_node=50, 
                            k_mnn=0, 
                            prune=False, 
                            device=device)
    
    info_log.print('--------> Integration graph is completed.')
    
    info_log.print('--------> Build KNN_adj for loss function ...')
    
    build_integration_loss_adj(adata, use_rep=use_rep, k=50, device=device)
    
    info_log.print('--------> KNN_adj is completed.')
    
    integration_gcn(adata, 
                    use_rep=use_rep,
                    max_epoch=max_epoch, 
                    lr=lr, 
                    num_features_gcn=32,
                    num_layers_gcn=2,  
                    device=device)

