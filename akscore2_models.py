from torch_geometric.nn import global_mean_pool
import torch
import torch.nn.functional as F
from torch import nn
from torch.nn import Linear
from torch_geometric.nn import GATv2Conv

class Akscore2_NonDock(torch.nn.Module):
    def __init__(self, node_dim=72, edge_dim=12, num_layers=5, hidden_dim=128, dropout=0.1):
        super().__init__()
        self.num_layers = num_layers
        self.node_dim = node_dim
        self.edge_dim = edge_dim
        self.hidden_dim = hidden_dim

        self.node_to_hidden = Linear(self.node_dim, self.hidden_dim)
        self.dropout_layer = nn.Dropout(dropout)

        self.protein_gnn_layers = []
        for num in range(self.num_layers):
            self.protein_gnn_layers.append(GATv2Conv(self.hidden_dim, self.hidden_dim, edge_dim=self.edge_dim))
        self.protein_gnn_layers = nn.ModuleList(self.protein_gnn_layers)

        self.ligand_gnn_layers = []
        for num in range(self.num_layers//2+1):
            self.ligand_gnn_layers.append(GATv2Conv(self.hidden_dim, self.hidden_dim, edge_dim=self.edge_dim))
        self.ligand_gnn_layers = nn.ModuleList(self.ligand_gnn_layers)

        self.graph_dec_rmsd = nn.Sequential(nn.Linear(self.hidden_dim*2, self.hidden_dim*2),
                                       nn.ReLU(),
                                       nn.Dropout(dropout),
                                       nn.Linear(self.hidden_dim*2, self.hidden_dim),
                                       nn.ReLU(),
                                       nn.Dropout(dropout),
                                       nn.Linear(self.hidden_dim, 1),)


    def forward(self, protein_graph_batch, ligand_graph_batch):

        protein_x, protein_edge_index, protein_edge_attr, protein_batch = protein_graph_batch.x, protein_graph_batch.edge_index, protein_graph_batch.edge_attr, protein_graph_batch.batch
        ligand_x, ligand_edge_index, ligand_edge_attr, ligand_batch = ligand_graph_batch.x, ligand_graph_batch.edge_index, ligand_graph_batch.edge_attr, ligand_graph_batch.batch

        protein_x = self.node_to_hidden(protein_x)
        ligand_x = self.node_to_hidden(ligand_x)

        for layer in self.protein_gnn_layers:
            protein_x = layer(protein_x, protein_edge_index, protein_edge_attr)
            protein_x = F.relu(protein_x)
            protein_x = self.dropout_layer(protein_x)

        for layer in self.ligand_gnn_layers:
            ligand_x = layer(ligand_x, ligand_edge_index, ligand_edge_attr)
            ligand_x = F.relu(ligand_x)
            ligand_x = self.dropout_layer(ligand_x)

        protein_x_feat = global_mean_pool(protein_x, protein_batch)  # [batch_size, hidden_channels]
        ligand_x_feat = global_mean_pool(ligand_x, ligand_batch)  # [batch_size, hidden_channels]
        protein_ligand_x_cat = torch.concat((protein_x_feat, ligand_x_feat), dim=-1)
        rmsd_logits = self.graph_dec_rmsd(protein_ligand_x_cat)

        return rmsd_logits




class Akscore2_DockS(torch.nn.Module):
    def __init__(self, node_dim=73, edge_dim=24, num_layers=5, hidden_dim=128, dropout=0.1):
        super().__init__()
        self.num_layers = num_layers
        self.node_dim = node_dim
        self.edge_dim = edge_dim
        self.hidden_dim = hidden_dim

        self.node_to_hidden = Linear(self.node_dim, self.hidden_dim)
        self.dropout_layer = nn.Dropout(dropout)

        self.protein_ligand_gnn_layers = []
        for num in range(self.num_layers):
            self.protein_ligand_gnn_layers.append(GATv2Conv(self.hidden_dim, self.hidden_dim, edge_dim=self.edge_dim))
        self.protein_ligand_gnn_layers = nn.ModuleList(self.protein_ligand_gnn_layers)

        self.graph_dec_bind = nn.Sequential(nn.Linear(self.hidden_dim, self.hidden_dim),
                                             nn.ReLU(),
                                             nn.Dropout(dropout),
                                             nn.Linear(self.hidden_dim, self.hidden_dim // 2),
                                             nn.ReLU(),
                                             nn.Dropout(dropout),
                                             nn.Linear(self.hidden_dim // 2, 1),
                                             )

        self.graph_dec_rmsd = nn.Sequential(nn.Linear(self.hidden_dim, self.hidden_dim),
                                             nn.ReLU(),
                                             nn.Dropout(dropout),
                                             nn.Linear(self.hidden_dim, self.hidden_dim // 2),
                                             nn.ReLU(),
                                             nn.Dropout(dropout),
                                             nn.Linear(self.hidden_dim // 2, 1),
                                             )

    def forward(self, graph_batch):

        x, edge_index, edge_attr, batch = graph_batch.x, graph_batch.edge_index, graph_batch.edge_attr, graph_batch.batch
        x = self.node_to_hidden(x)

        for layer in self.protein_ligand_gnn_layers:
            x = layer(x, edge_index, edge_attr)
            x = F.relu(x)
            x = self.dropout_layer(x)

        protein_ligand_x_dock = global_mean_pool(x, batch)  # [batch_size, hidden_channels]
        bind_logits = self.graph_dec_bind(protein_ligand_x_dock)
        rmsd_logits = self.graph_dec_rmsd(protein_ligand_x_dock)

        return bind_logits, rmsd_logits


class Akscore2_DockC(torch.nn.Module):
    def __init__(self, node_dim=73, edge_dim=24, num_layers=5, hidden_dim=256, dropout=0.1):
        super().__init__()
        self.num_layers = num_layers
        self.node_dim = node_dim
        self.edge_dim = edge_dim
        self.hidden_dim = hidden_dim

        self.node_to_hidden = Linear(self.node_dim, self.hidden_dim)
        self.dropout_layer = nn.Dropout(dropout)

        self.protein_ligand_gnn_layers = []
        for num in range(self.num_layers):
            self.protein_ligand_gnn_layers.append(GATv2Conv(self.hidden_dim, self.hidden_dim, edge_dim=self.edge_dim))
        self.protein_ligand_gnn_layers = nn.ModuleList(self.protein_ligand_gnn_layers)

        self.graph_dec_bind = nn.Sequential(nn.Linear(self.hidden_dim, self.hidden_dim),
                                            nn.ReLU(),
                                            nn.Dropout(dropout),
                                            nn.Linear(self.hidden_dim, self.hidden_dim // 2),
                                            nn.ReLU(),
                                            nn.Dropout(dropout),
                                            nn.Linear(self.hidden_dim // 2, 1),
                                            )

    def forward(self, graph_batch):

        x, edge_index, edge_attr, batch = graph_batch.x, graph_batch.edge_index, graph_batch.edge_attr, graph_batch.batch
        x = self.node_to_hidden(x)

        for layer in self.protein_ligand_gnn_layers:
            x = layer(x, edge_index, edge_attr)
            x = F.relu(x)
            x = self.dropout_layer(x)

        protein_ligand_x_dock = global_mean_pool(x, batch)  # [batch_size, hidden_channels]
        bind_logits = self.graph_dec_bind(protein_ligand_x_dock)

        return bind_logits

