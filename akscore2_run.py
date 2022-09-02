import os, pickle
import numpy as np
import argparse

from sklearn.metrics import mean_squared_error, mean_absolute_error, auc
from torch_geometric.nn import global_mean_pool
from torch_geometric.loader import DataLoader
from torch_geometric.data import Data

import torch
import torch.nn.functional as F
from torch import nn
from torch.nn import Linear
from torch_geometric.nn import GATv2Conv
from torch_geometric.data import Dataset
from torch_geometric.data import batch as pyg_batch_func
from sklearn.metrics import roc_auc_score

import random
import time
import gzip
from itertools import chain

from construct_pyg_graph import pdb_list_cut, make_graph
from akscore2_models import Akscore2_DockC, Akscore2_DockS, Akscore2_NonDock


class akscore2_dataset(Dataset):
    def __init__(self, receptor_path, ligand_path):
        super(akscore2_dataset, self).__init__()
        self.receptor_path =receptor_path
        self.ligand_txt_list = pdb_list_cut(ligand_path)
    def len(self):
        return len(self.ligand_txt_list)

    def graph_modification(self, graph):
        x, edge_index, edge_attr = graph.x.detach().clone(), graph.edge_index.detach().clone(), graph.edge_attr.detach().clone()

        protein_edge_attr_idx = torch.where((edge_attr[:, :3] == torch.Tensor([1, 0, 0])).all(dim=1))[0]
        ligand_edge_attr_idx = torch.where((edge_attr[:, :3] == torch.Tensor([0, 1, 0])).all(dim=1))[0]

        #### remove docking feature for platform app
        x = torch.concat((x[:, :-5], x[:, -4:]), axis=1)
        edge_attr = edge_attr[:, 3:-9]

        protein_edge_index = edge_index[:, protein_edge_attr_idx]
        ligand_edge_index = edge_index[:, ligand_edge_attr_idx]

        protein_ligand_node_sep_idx = torch.min(ligand_edge_index)

        protein_x = x[:protein_ligand_node_sep_idx, :]
        ligand_x = x[protein_ligand_node_sep_idx:, :]

        protein_graph = Data(x=protein_x, edge_index=protein_edge_index, edge_attr=edge_attr[protein_edge_attr_idx, :])
        ligand_graph = Data(x=ligand_x, edge_index=ligand_edge_index - torch.min(ligand_edge_index), edge_attr=edge_attr[ligand_edge_attr_idx, :])



        return graph, protein_graph, ligand_graph

    def get(self, idx):

        graph, error_graph_tag = make_graph(self.receptor_path, self.ligand_txt_list[idx])
        graph, protein_graph, ligand_graph = self.graph_modification(graph)

        return graph, protein_graph, ligand_graph


if __name__ == "__main__" :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-r','--receptor_path', default='input_example/1a30_protein.pdb', help='receptor pdb')
    parser.add_argument('-l','--ligand_path', default='input_example/test_pdb.pdb', help='ligands dlg or pdbqt or dlg, pdbqt list')
    parser.add_argument('-o','--output', default='result.csv', help='result output file')
    parser.add_argument('-s','--select_dock_model', default='akscore2_dockc', help='select either akscore2_dockc or akscore_docks')


    parser.add_argument('-ndw','--akscore2_nondock_weight_path', default='model_weights/akscore2_nondock_weights.pth', help='akscore2_nondock_weight_path')
    parser.add_argument('-dsw','--akscore2_docks_weight_path', default='model_weights/akscore2_docks_weights.pth', help='akscore2_docks_weight_path')
    parser.add_argument('-dcw','--akscore2_dockc_weight_path', default='model_weights/akscore2_dockc_weights.pth', help='akscore2_dockc_weight_path')

    parser.add_argument('--batch_size', default=8, type=int, help='batch size')
    parser.add_argument('--num_workers', default=0, type=int, help="cpu worker number")
    parser.add_argument('--device', type=str, default='gpu', help='choose device: cpu or gpu')

    args = parser.parse_args()


    if args.device == 'cpu':
        device = torch.device("cpu")
    else:
        if torch.cuda.is_available():
            device = torch.device("cuda:0")
        else:
            print("gpu is not available, run on cpu")
            device = torch.device("cpu")

    model_nondock = Akscore2_NonDock().to(device)
    ch_nondock = torch.load(args.akscore2_nondock_weight_path, map_location=device)

    if args.select_dock_model == "akscore2_dockc":
        model_dock = Akscore2_DockC().to(device)
        ch_dock = torch.load(args.akscore2_dockc_weight_path, map_location=device)

    else:
        model_dock = Akscore2_DockS().to(device)
        ch_dock = torch.load(args.akscore2_docks_weight_path, map_location=device)

    model_nondock.load_state_dict(ch_nondock['model'])
    model_dock.load_state_dict(ch_dock['model'])

    model_nondock.eval().to(device)
    model_dock.eval().to(device)


    dataset = akscore2_dataset(args.receptor_path, args.ligand_path)
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=False, num_workers=args.num_workers)

    with torch.no_grad():
        for idx, (graph_batch, protein_graph_batch, ligand_graph_batch) in enumerate(loader):
            graph_batch.to(device)
            protein_graph_batch.to(device)
            ligand_graph_batch.to(device)
            print(f"{idx}/{len(loader)}")
            pred_nondock = model_nondock(protein_graph_batch, ligand_graph_batch)

            if args.select_dock_model == "akscore2_dockc":
                pred_dock_b = model_dock(graph_batch)
            else:
                pred_dock_b, pred_dock_r = model_dock(graph_batch)
                pred_dock_b = pred_dock_b+pred_dock_r

            pred_nondock_sigmoid = torch.sigmoid(pred_nondock)
            pred_ens = pred_nondock_sigmoid*pred_dock_b
            pred_ens = pred_ens.cpu().detach().numpy()
            print(pred_ens)
