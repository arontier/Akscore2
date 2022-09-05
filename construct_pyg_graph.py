#!/usr/bin/env python


import sys, os, pickle, gzip
import numpy as np
import pandas as pd
import subprocess

from tqdm import tqdm
from glob import glob
from multiprocessing import Pool

import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
import rdkit.Chem.rdchem as rdchem
from rdkit.Chem.rdchem import BondStereo, BondType, HybridizationType
from rdkit.Chem import rdmolfiles
from rdkit.Chem import rdPartialCharges
from rdkit import RDLogger 
#RDLogger.DisableLog('rdApp.*')

import torch
from torch_geometric.data import Data
#from torch_geometric.loader import DataLoader

import logging
from meeko import PDBQTMolecule, MoleculePreparation, rdkitutils
from openbabel import pybel as pb
from io import StringIO, BytesIO
from rdkit.Chem.rdmolfiles import ForwardSDMolSupplier
import random
import math
import random
from time import perf_counter
from openbabel.openbabel import *



logger = logging.getLogger()
logger.setLevel(logging.INFO)

#node_attr = Symbol[20] + Degeree[7] + NumHs[5] + ImpVal[6] + ExpVal[6] + ExpHs[2] + ImpHs[5] + Hybri[5] + Arom[1] + IsRing + Radical[3] + FormalCharge[7] + Hydropho + donor + acceptor + acidic + basic = [73]
#edge_attr = PLbond[3] + bondtype[4] + Streo[6] + Ring + conjugate + PLbondtype[5] + bond_distance[4]= [24]
logging.basicConfig(level=logging.INFO)

os.environ['CUDA_LAUNCH_BLOCKING'] = "1"
os.environ["CUDA_VISIBLE_DEVICES"] = "0"


prot_atom_hydrophobicity = {
 'MET': {'C':1, 'CA':0, 'CB':0, 'CE':0, 'CG':0, 'H':1, 'HA':0, 'HB1':0, 'HB2':0, 'HE1':0, 'HE2':0, 'HE3':0, 'HG1':0, 'HG2':0, 'N':1, 'O':1, 'SD':1, 'OXT':1, 'HN1':1, 'HN2':1, 'HN3':1},
 'ASP': {'C':1,  'CA':0,  'CB':0,  'CG':1,  'H':1,  'HA':0,  'HB1':0,  'HB2':0,  'N':1,  'O':1,  'OD1':1,  'OD2':1, 'OXT':1, 'HN1':1, 'HN2':1, 'HN3':1},
 'PRO': {'C':1,  'CA':0,  'CB':0,  'CD':0,  'CG':0,  'HA':0,  'HB1':0,  'HB2':0,  'HD1':0,  'HD2':0,  'HG1':0,  'HG2':0,  'N':0,  'O':1, 'OXT':1, 'HN1':1, 'HN2':1, 'HN3':1},
 'TYR': {'C':1,  'CA':0,  'CB':0,  'CD1':0,  'CD2':0,  'CE1':0,  'CE2':0,  'CG':0,  'CZ':0,  'H':1,  'HA':0,  'HB1':0,  'HB2':0,  'HD1':0,  'HD2':0,  'HE1':0,  'HE2':0,  'HH':1,  'N':1,  'O':1,  'OH':1, 'OXT':1, 'HN1':1, 'HN2':1, 'HN3':1},
 'LEU': {'1HD1':0,  '1HD2':0,  '2HD1':0,  '2HD2':0,  '3HD1':0,  '3HD2':0,  'C':1,  'CA':0,  'CB':0,  'CD1':0,  'CD2':0,  'CG':0,  'H':1,  'HA':0,  'HB1':0,  'HB2':0,  'HG':0,  'N':1,  'O':1, 'OXT':1, 'HN1':1, 'HN2':1, 'HN3':1},
 'VAL': {'1HG1':0,  '1HG2':0,  '2HG1':0,  '2HG2':0,  '3HG1':0,  '3HG2':0,  'C':1,  'CA':0,  'CB':0,  'CG1':0,  'CG2':0,  'H':1,  'HA':0,  'HB':0,  'N':1,  'O':1, 'OXT':1, 'HN1':1,'HN2':1,'HN3':1},
 'GLY': {'C':1, 'CA':0, 'H':1, 'HA1':0, 'HA2':0, 'N':1, 'O':1, 'OXT':1,'HN1':1,'HN2':1,'HN3':1},
 'ASN': {'1HD2':1,  '2HD2':1,  'C':1,  'CA':0,  'CB':0,  'CG':1,  'H':1,  'HA':0,  'HB1':0,  'HB2':0,  'HN1':1,  'HN2':1,  'HN3':1,  'N':1,  'ND2':1,  'O':1,  'OD1':1, 'OXT':1},
 'THR': {'1HG2':0, '2HG2':0,  '3HG2':0,  'C':1,  'CA':0,  'CB':0,  'CG2':0,  'H':1,  'HA':0,  'HB':0,  'HG1':1,  'N':1,  'O':1,  'OG1':1, 'OXT':1, 'HN1':1, 'HN2':1, 'HN3':1},
 'HIS': {'C':1,  'CA':0,  'CB':0,  'CD2':0,  'CE1':1,  'CG':0,  'H':1,  'HA':0,  'HB1':0,  'HB2':0,  'HD2':0,  'HE1':1,  'HN1':1,  'HN2':1,  'HN3':1,  'N':1,  'ND1':1,  'NE2':1,  'O':1, 'OXT':1},
 'SER': {'C':1, 'CA':0, 'CB':0, 'H':1, 'HA':0, 'HB1':0, 'HB2':0, 'HG':1, 'N':1, 'O':1, 'OG':1, 'OXT':1, 'HN1':1, 'HN2':1, 'HN3':1},
 'PHE': {'C':1,  'CA':0,  'CB':0,  'CD1':0,  'CD2':0,  'CE1':0,  'CE2':0,  'CG':0,  'CZ':0,  'H':1,  'HA':0,  'HB1':0,  'HB2':0,  'HD1':0,  'HD2':0,  'HE1':0,  'HE2':0,  'HN1':1,  'HN2':1,  'HN3':1,  'HZ':0, 'N':1,  'O':1, 'OXT':1},
 'ILE': {'1HD1':0,  '1HG1':0,  '1HG2':0,  '2HD1':0,  '2HG1':0,  '2HG2':0,  '3HD1':0,  '3HG2':0,  'C':1,  'CA':0,  'CB':0,  'CD1':0,  'CG1':0,  'CG2':0,  'H':1,  'HA':0,  'HB':0,  'N':1,  'O':1, 'OXT':1,'HN1':1,'HN2':1,'HN3':1},
 'GLU': {'C':1,  'CA':0,  'CB':0,  'CD':1,  'CG':0,  'H':1,  'HA':0,  'HB1':0,  'HB2':0,  'HG1':0,  'HG2':0,  'N':1,  'O':1,  'OE1':1,  'OE2':1, 'OXT':1,'HN1':1,'HN2':1,'HN3':1},
 'LYS': {'C':1,  'CA':0,  'CB':0,  'CD':0,  'CE':0,  'CG':0,  'H':1,  'HA':0,  'HB1':0,  'HB2':0,  'HD1':0,  'HD2':0,  'HE1':0,  'HE2':0,  'HG1':0,  'HG2':0,  'HZ1':1,  'HZ2':1,  'HZ3':1,  'N':1,  'NZ':1,  'O':1, 'OXT':1, 'HN1':1, 'HN2':1, 'HN3':1},
 'ARG': {'1HH1':1,  '1HH2':1,  '2HH1':1,  '2HH2':1,  'C':1,  'CA':0,  'CB':0,  'CD':0,  'CG':0,  'CZ':1,  'H':1,  'HA':0,  'HB1':0,  'HB2':0,  'HD1':0,  'HD2':0,  'HE':1,  'HG1':0,  'HG2':0, 'HN1':1,  'HN2':1,  'HN3':1,  'N':1,  'NE':1,  'NH1':1,  'NH2':1,  'O':1, 'OXT':1},
 'TRP': {'C':1,  'CA':0,  'CB':0,  'CD1':0,  'CD2':0,  'CE2':0,  'CE3':0,  'CG':0,  'CH2':0,  'CZ2':0,  'CZ3':0,  'H':1,  'HA':0,  'HB1':0,  'HB2':0,  'HD1':0,  'HE1':1,  'HE3':0,  'HH2':1,  'HZ2':0,  'HZ3':0,  'N':1,  'NE1':1,  'O':1, 'OXT':1, 'HN1':1, 'HN2':1, 'HN3':1},
 'GLN': {'1HE2':1,  '2HE2':1,  'C':1,  'CA':0,  'CB':0,  'CD':1,  'CG':0,  'H':1,  'HA':0,  'HB1':0,  'HB2':0,  'HG1':0,  'HG2':0,  'HN1':1,  'HN2':1,  'HN3':1,  'N':1,  'NE2':1,  'O':1,  'OE1':1, 'OXT':1},
 'ALA': {'C':1, 'CA':0, 'CB':0, 'H':1, 'HA':0, 'HB1':0, 'HB2':0, 'HB3':0, 'N':1, 'O':1, 'OXT':1, 'HN1':1, 'HN2':1, 'HN3':1},
 'CYS': {'C':1, 'CA':0, 'CB':0, 'H':1, 'HA':0, 'HB1':0, 'HB2':0, 'HG':1, 'N':1, 'O':1, 'SG':1, 'OXT':1, 'HN1':1, 'HN2':1, 'HN3':1},
 'CSD': {'C':1, 'CA':0, 'CB':0, 'H':1,  'HA':0,  'HB1':0,  'HB2':0,  'N':1,  'O':1,  'OD1':1,  'OD2':1,  'SG':1, 'OXT':1},
 'MLY': {'C':1, 'CA':0, 'CB':0, 'CD':0, 'CE':0, 'CG':0, 'CH1':0, 'CH2':0, 'H':1, 'HZ':0, 'N':1, 'NZ':1, 'O':1, 'OXT':1},
 'PCA': {'C':1, 'CA':0, 'CB':0, 'CD':0, 'CG':0, 'N':1, 'O':1, 'OE':1, 'OXT':1},
 'PTR': {'C':1, 'CA':0, 'CB':0, 'CD1':0, 'CD2':0, 'CE1':0, 'CE2':0, 'CG':0, 'CZ':0, 'H':1, 'N':1, 'O':1, 'O1P':1,'O2P':1, 'O3P':1, 'OH':1, 'P':1, 'OXT':1},
 'CSO': {'C':1, 'CA':0, 'CB':0, 'H':1, 'HA':0, 'HB1':0, 'HB2':0, 'HD':0, 'N':1, 'O':1, 'OD':1, 'SG':1, 'OXT':1},
 'KCX': {'C':1, 'CA':0, 'CB':0, 'CD':0, 'CE':0, 'CG':0, 'CX':0, 'H':1, 'HZ':1, 'N':1, 'NZ':1, 'O':1, 'OQ1':1, 'OQ2':1, 'OXT':1},
 'ACE': {'C':1, 'CA':0, 'HA1':0, 'HA2':0, 'HA3':0, 'O':1, 'OXT':1},
 'MSE': {'C':1, 'CA':0, 'CB':0, 'CE':0, 'CG':0, 'H':1, 'N':1, 'O':1, 'SE':1, 'OXT':1,'HN1':1,'HN2':1,'HN3':1},
 'TPO': {'C':1, 'CA':0, 'CB':0, 'CG2':0, 'H':1, 'N':1, 'O':1, 'O1P':1, 'O2P':1, 'O3P':1, 'OG1':1, 'P':1, 'OXT':1},
 'SEP': {'C':1, 'CA':0, 'CB':0, 'H':1, 'N':1, 'O':1, 'O1P':1, 'O2P':1, 'O3P':1, 'OG':1, 'P':1, 'OXT':1},
 'LLP': {'C':1, "C2'":0, 'C2':0, 'C3':1, "C4'":0, 'C4':0, "C5'":0, 'C5':0, 'C6':0, 'CA':0, 'CB':0, 'CD':0, 'CE':0, 'CG':0, 'H':1, 'H21':0, 'H22':0, 'H23':0, 'H41':0, 'H42':0, 'H51':0, 'H52':0, 'H6':0, 'HA':0, 'HB1':0, 'HB2':0, 'HD1':0, 'HD2':0, 'HE1':0, 'HE2':0, 'HG1':0, 'HG2':0, 'HZ1':1, 'HZ2':1, 'N':1, 'N1':0, 'NZ':1, 'O':1, 'O3':1,  'OP1':1, 'OP2':1, 'OP3':1, 'OP4':1, 'P':1, 'OXT':1},
 'SO4': {'O1':1, 'O2':1, 'O3':1, 'O4':1, 'S':1, 'OXT':1},
 'HOH': {'O':1},
 'NA': {'NA':1},
 'ZN': {'ZN':1},
 'MG': {'MG':1},
 'K': {'K':1},
 'MN': {'MN':1},
 'CA': {'CA':1},
 'CD': {'CD':1},
 'NI': {'NI':1},
 'CU': {'CU':1},
 'CAS': {'CA':1},
 'HG': {'HG':1},
 'CO': {'CO':1},
 'FE': {'FE':1},
 'FE2': {'FE':1},
 'CS': {'CS':1},
 'YCM': {'OXT':1},
 'LI': {'LI':1},
 'IAS': {'OXT':1},
 'SR': {'SR':1},
 'CMT': {'OXT':1},
 'DAL': {'OXT':1},
 'CAF': {'CA':1},
 'CU1': {'CU':1},
 'GA': {'GA':1},
 'RB': {'RB':1},
 'AU': {'AU':1},
 'IN': {'IN':1},

}

def one_hot(x, allowable_set):
    if x not in allowable_set:            
        x = allowable_set[-1]
    return list(map(lambda s: x == s, allowable_set))

def Mol22mol(file=None):
    mols=[]
    obconversion = OBConversion()
    obconversion.SetInAndOutFormats("mol2","pdb")
    obmol = OBMol()

    notatend = obconversion.ReadFile(obmol,file)
    while notatend:
        mols.append(Chem.MolFromPDBBlock(obconversion.WriteString(obmol)))
        obmol = OBMol()
        notatend = obconversion.Read(obmol)
    return mols


def pdb2mol(pdb):
    try:
        return Chem.MolFromPDBFile(pdb)
    except:
        return Chem.MolFromPDBBlock(pdb,sanitize=True)


def pdb2data(file):
    pocket_mol = pdb2mol(file)
    pocket_heavy    = []
    pocket_hydrogen = []    
    bond = []    
    xyz_pdb = {}
    try:
        for line in open(file).readlines():
            if line[0:4] in ['ATOM'] and line.split()[3] == 'H':
                pocket_hydrogen.append(line)            
            elif line[0:4] in ['ATOM'] and line.split()[-1] != 'H':
                pocket_heavy.append(line)
            elif line[0:4] in ['HETA'] and line.split()[3] != 'HOH':
                pocket_heavy.append(line)            
            elif line[0:4] in ['CONE']:
                bond.append( line.split() )  
    except:
        for line in file.split("\n"):
            if line[0:4] in ['ATOM'] and line.split()[3] == 'H':
                pocket_hydrogen.append(line)            
            elif line[0:4] in ['ATOM'] and line.split()[-1] != 'H':
                pocket_heavy.append(line)
            elif line[0:4] in ['HETA'] and line.split()[3] != 'HOH':
                pocket_heavy.append(line)            
            elif line[0:4] in ['CONE']:
                bond.append( line.split() )  
    
    name = file[:-4]


    pocket_heavy = { idx: [ line[12:16].strip(), line[17:20].strip(), [ float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()) ] ,idx] for idx, line in enumerate(pocket_heavy) }
    pocket_hydro = { int(line[6:11]): [ line[12:16].strip(), line[17:20].strip(), [ float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip()) ] ] for idx, line in enumerate(pocket_hydrogen) } 
    bond   = { int(line[1]): [ int(num) for num in line[2:] ] for line in bond }
        
    return pocket_mol, pocket_heavy, bond


def get_default_config():
    config = {
            'keep_nonpolar_hydrogens': False,
            'hydrate': False,
            'flexible_amides': True,
            'rigid_macrocycles': False,
            'min_ring_size': 7,
            'max_ring_size': 33,
            'rigidify_bonds_smarts': [],
            'rigidify_bonds_indices': [],
            'double_bond_penalty': 50,
            'atom_type_smarts': {},
            'add_index_map': False,
            'remove_smiles': False,
            }
    return config


def rdkit_mol_TO_pdbqt(mol, config=None) :
    if not config :
        config = get_default_config()

    mol_counter = 0
    num_skipped = 0

    is_valid = mol is not None
    mol_counter += int(is_valid==True)
    num_skipped += int(is_valid==False)
    if not is_valid:
        return

    preparator = MoleculePreparation.from_config(config)
    preparator.prepare(mol)

    try :
        ligand_prepared = preparator.write_pdbqt_string()
    except RuntimeError :
        return

    return ligand_prepared


def poss2data(poss_mols,clus_list=None):
    datas = []
    lig_mol = []
    dock_score = []
    for i_mol,mol in enumerate(poss_mols) :
        try :
            pdbqt_string=mol.write_pdbqt_string()
            mol = mol.export_rdkit_mol(exclude_hydrogens=False) # 수소 위치가 필요 없으면 True
        except RuntimeError :
            pdbqt_string=mol.write_pdbqt_string()
            obmol = pb.readstring("pdbqt",pdbqt_string)
            sdf_string=obmol.write("sdf")
    
            str_io = StringIO(sdf_string)
            str_data = str_io.read().encode('utf8')
            bytes_io = BytesIO(str_data)
            suppl = ForwardSDMolSupplier(bytes_io,removeHs=False)
            for mol in suppl:
                break
    
        Chem.AssignStereochemistryFrom3D(mol)
#        mol = Chem.RemoveHs(mol)
        l = 0
        data = {}
        xyz_pdbqt = {}
        #from pdbqt
        for line in pdbqt_string.split("\n"):
            if line.startswith("ATOM"):
                xyz_pdbqt[float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())] = float(line[71:77].strip()),line[77:79].strip()
        
        #from mol
        conf = mol.GetConformer()
        nums = 0
        for i, atom in enumerate(mol.GetAtoms()):
            atom_xyz = list(conf.GetAtomPosition(i))
            if (atom_xyz[0],atom_xyz[1],atom_xyz[2]) in xyz_pdbqt:
                data[i] = [atom.GetSymbol(),xyz_pdbqt[atom_xyz[0],atom_xyz[1],atom_xyz[2]][0], atom_xyz,xyz_pdbqt[atom_xyz[0],atom_xyz[1],atom_xyz[2]][1],nums]
                nums += 1

    return mol, data


def calculateDistance(arr1, arr2):
    return np.linalg.norm( arr1[:, None, :] - arr2[None, :, :], axis = -1)


class Featurize:
    def __init__(self):
        self.H_donor    = Chem.MolFromSmarts("[$([N;!H0;v3,v4&+1]),$([O,S;H1;+0]),n&H1&+0]")
        self.H_acceptor = Chem.MolFromSmarts("[$([O,S;H1;v2;!$(*-*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=[O,N,P,S])]),n&H0&+0,$([o,s;+0;!$([o,s]:n);!$([o,s]:c:n)])]")
        self.acidic     = Chem.MolFromSmarts("[$([C,S](=[O,S,P])-[O;H1,-1])]")
        self.basic      = Chem.MolFromSmarts("[#7;+,$([N;H2&+0][$([C,a]);!$([C,a](=O))]),$([N;H1&+0]([$([C,a]);!$([C,a](=O))])[$([C,a]);!$([C,a](=O))]),$([N;H0&+0]([C;!$(C(=O))])([C;!$(C(=O))])[C;!$(C(=O))])]")        
    
    def AtomNode(self, mol, data, pocket=True):
        # AllChem.ComputeGasteigerCharges(mol)
        H_donor_match    = sum(mol.GetSubstructMatches(self.H_donor), ())
        H_acceptor_match = sum(mol.GetSubstructMatches(self.H_acceptor), ())
        acidic_match     = sum(mol.GetSubstructMatches(self.acidic), ())
        basic_match      = sum(mol.GetSubstructMatches(self.basic), ())
        
        node = []
        xyz  = []                   
        for idx, atom in enumerate(mol.GetAtoms()): 
            if idx in data:
                if pocket:
                    try:
                        hydrophobicity = prot_atom_hydrophobicity[ data[idx][1] ][ data[idx][0] ]
                    except:
                        hydrophobicity = 0 
                else:
                    hydrophobicity = 0 if abs(data[idx][1]) <0.2 else 1
                    
                tmp = np.array(
                    one_hot( atom.GetSymbol(), ['C', 'H', 'O', 'N', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'B', 'Se', 'Fe', 'Ru', 'Mn', 'Co', 'Ni', 'Cu', 'Zn', 'other']) +
                    one_hot( atom.GetDegree(), [0, 1, 2, 3, 4, 5, 6]) +                
                    one_hot( atom.GetTotalNumHs(), [0, 1, 2, 3, 4]) +                
                    one_hot( atom.GetImplicitValence(), [0, 1, 2, 3, 4, 5]) +                
                    one_hot( atom.GetExplicitValence(), [0, 1, 2, 3, 4, 5]) +                
                    one_hot( atom.GetNumExplicitHs(), [0, 1]) +                
                    one_hot( atom.GetNumImplicitHs(), [0, 1, 2, 3, 4]) +                
                    one_hot( atom.GetHybridization(), [ HybridizationType.SP, HybridizationType.SP2, HybridizationType.SP3, HybridizationType.SP3D, HybridizationType.SP3D2] ) +
                    one_hot( atom.GetFormalCharge(), [-2,-1,0,1,2,3,4]) +
                    one_hot( atom.GetNumRadicalElectrons(), [0,1,2]) +

                    [ atom.GetIsAromatic(), atom.IsInRing(), hydrophobicity ] +
                    [ idx in H_donor_match, idx in H_acceptor_match, idx in acidic_match, idx in basic_match ]            
                )
                node.append(tmp)                       
                xyz.append(data[ idx ][2])
            
        return node, xyz
    
        
    def CommonEdge(self, pocket=False, ligand=False):
        if pocket == True and ligand == False: ## pocket covalent bond
            return [1, 0, 0] 
        elif pocket == False and ligand == True: ## ligand covalent bond
            return [0, 1, 0] 
        elif pocket == True and ligand == True: ## PL non-covalent bond
            return [0, 0, 1]        
    
    def CovalentEdge(self, mol,data, pocket=False, ligand=False, num=0):
        commonEdge = self.CommonEdge(pocket=pocket, ligand=ligand)            
        edge_index, edge_attr = [], []  
        for idx, bond in enumerate(mol.GetBonds()):
            if bond.GetBeginAtomIdx() in data and bond.GetEndAtomIdx() in data:
                edge_index.append( [ data[bond.GetBeginAtomIdx()][-1] + num, data[bond.GetEndAtomIdx()][-1] + num ] )
                edge_index.append( [ data[bond.GetEndAtomIdx()][-1] + num, data[bond.GetBeginAtomIdx()][-1] + num ] )            
                tmp = np.array( 
                    commonEdge +
                    one_hot( bond.GetBondType(), [ BondType.SINGLE, BondType.DOUBLE, BondType.TRIPLE, rdchem.BondType.AROMATIC ] ) +
                    one_hot( bond.GetStereo(),   [ BondStereo.STEREOANY, BondStereo.STEREOCIS, BondStereo.STEREOE, BondStereo.STEREONONE, BondStereo.STEREOTRANS, BondStereo.STEREOZ ] ) +
                    [ bond.IsInRing(), bond.GetIsConjugated() ] +         
                    [0, 0, 0, 0, 0] + [1, 1, 1, 1]
                )
                edge_attr.append(tmp)
                edge_attr.append(tmp)   
            
        return edge_index, edge_attr
    
    def PL_NonCovalentEdge(self, pocket_node, ligand_node, distance_matrix):
        num = len(pocket_node)
        distance_4 = np.array( np.where( distance_matrix <= 8 ) )
        commonEdge = self.CommonEdge(pocket=True, ligand=True)
        
        edge_index, edge_attr = [], []
        for idx, prot_idx in enumerate(distance_4[0]):
            
            lig_idx = distance_4[1][idx]
            lig_idx_total = lig_idx + num
        
            edge_index.append( [ prot_idx, lig_idx_total ] )
            edge_index.append( [ lig_idx_total, prot_idx ] )            
       
            Hydrophobic = pocket_node[prot_idx][-5], ligand_node[lig_idx][-5] ## if (0,0) == hydrophobic bond
            Hydrophobic_bond = 1 if Hydrophobic == (0, 0) else 0                
            
            PL_Hydrogen = pocket_node[prot_idx][-4], ligand_node[lig_idx][-3]
            PL_Hydrogen_bond = 1 if PL_Hydrogen == (1, 1) else 0
            
            LP_Hydrogen = pocket_node[prot_idx][-3], ligand_node[lig_idx][-4]
            LP_Hydrogen_bond = 1 if LP_Hydrogen == (1, 1) else 0
            
            prot_ionic_prop = pocket_node[prot_idx][-6]
            lig_ionic_prop  = ligand_node[lig_idx][-6]
            
            Ionic_bond_PL = 1 if prot_ionic_prop < 0 and lig_ionic_prop > 0 else 0
            Ionic_bond_LP = 1 if prot_ionic_prop > 0 and lig_ionic_prop < 0 else 0

            if distance_matrix[prot_idx][lig_idx] < 2:
                ang_2,ang_4,ang_6,ang_8 = 1,1,1,1
            elif distance_matrix[prot_idx][lig_idx] < 4:
                ang_2,ang_4,ang_6,ang_8 = 0,1,1,1
            elif distance_matrix[prot_idx][lig_idx] < 6:
                ang_2,ang_4,ang_6,ang_8 = 0,0,1,1
            elif distance_matrix[prot_idx][lig_idx] < 8:
                ang_2,ang_4,ang_6,ang_8 = 0,0,0,1
            
            edge_tmp = commonEdge + [0] * 12 + [Hydrophobic_bond, PL_Hydrogen_bond, LP_Hydrogen_bond, Ionic_bond_PL, Ionic_bond_LP ] + [ang_2,ang_4,ang_6,ang_8]
            edge_attr.append(edge_tmp)
            edge_attr.append(edge_tmp)
            
        return edge_index, edge_attr
    

def mol_2_graph(ligand_pdb, protein_file,clus_list=None):
    featurize = Featurize()
    
    ## pocket file and ligand file to block data
    ## make rdkit.mol for node and edge
    try:
        ligand_mol = pdb2mol(ligand_pdb)
        ligand_poss = rdkit_mol_TO_pdbqt(ligand_mol)
        poss_mols = PDBQTMolecule(ligand_poss)

        ligand_mol, ligand_data = poss2data(poss_mols)
        ## node part
        if ligand_mol == None:
            raise Exception("no have ligand")
        pocket_file = make_pocket(protein_file,ligand_data)
        pocket_mol, pocket_data, bond = pdb2data(pocket_file)

        if pocket_mol == None:
            raise Exception("no have pocket")
        node_pocket, xyz_pocket = featurize.AtomNode(pocket_mol, pocket_data, pocket=True) 
        node_ligand, xyz_ligand = featurize.AtomNode(ligand_mol, ligand_data, pocket=False)
        node_attr = node_pocket + node_ligand
        node_attr = np.array(node_attr)
        P_num = len(node_pocket)
           
            
        ## calculate distance for Pocket-Ligand non-covalent bond        
        distance_matrix_PL = calculateDistance( np.array(xyz_pocket), np.array(xyz_ligand) )
        
        ## covalent edge ( Pocket-Pocket, Ligand-Ligand atoms)
        edge_index_pocket, edge_attr_pocket = featurize.CovalentEdge(pocket_mol, pocket_data, pocket=True, ligand=False, num=0)  # num == match index
        edge_index_ligand, edge_attr_ligand = featurize.CovalentEdge(ligand_mol, ligand_data, pocket=False, ligand=True, num=P_num)        
                
        ## non covalent edge ( Pocket - Ligand )
        edge_index_PL_non, edge_attr_PL_non = featurize.PL_NonCovalentEdge(node_pocket, node_ligand, distance_matrix_PL)
        
        ## total edge for Graph
        edge_index = edge_index_pocket + edge_index_ligand + edge_index_PL_non
        edge_index = np.array(edge_index)
        edge_attr  = edge_attr_pocket  + edge_attr_ligand  + edge_attr_PL_non
        edge_attr = np.array(edge_attr)
            
        node_attrs = torch.tensor(node_attr,dtype=torch.float)
        edge_attrs = torch.tensor(edge_attr,dtype=torch.float)
        edge_index = torch.tensor(edge_index,dtype=torch.long)
        edge_indexs = edge_index.t().contiguous()
        error = 1 
    except Exception as E:
        logging.info(f"{E}, error ligand")
        node_attrs = torch.tensor(np.zeros([100,73]),dtype=torch.float)
        edge_attrs = torch.tensor(np.zeros([1000,24]),dtype=torch.float)
        edge_index = torch.tensor(np.zeros([1000,2]),dtype=torch.long)
        edge_indexs = edge_index.t().contiguous()
        error = 0

    ## data to tensor for Graph
        
    return node_attrs, edge_attrs, edge_indexs, error


def cluster_list(dlg):
    RMSD_TABLE = '_____|______|______|___________|_________|______'
    line_s = []
    Energy_s = []
    with open(dlg) as in_dlg :
        rmsd_table_gate = False
        energy_table_gate = False
        for line in in_dlg:
            if line.startswith(RMSD_TABLE) :
                rmsd_table_gate = True
                continue
            if rmsd_table_gate and line.strip() == '':
                break
            if rmsd_table_gate:
                line_s.append(line)
            if line.startswith('DOCKED: MODEL') :
                energy_table_gate = True
                energy = []
            if line.startswith('DOCKED: ENDMDL') :
                energy_table_gate = False 
                Energy_s.append(energy[0])
            if energy_table_gate and line.startswith('DOCKED: USER'):
                if 'Energy' in line:
                    energy.append(float(line.split("=")[1].split()[0]))


    cluster_zero = [int(line.split()[2]) for line in line_s if line.split()[4].strip()=="0.00"]
    if len(cluster_zero) > 100:
        cluster_zero = cluster_zero[:100]
    return cluster_zero, Energy_s

def PDBQT_native(dlg):
    energy = []
    with open(dlg) as in_dlg :
        for line in in_dlg:
            if line.startswith("INPUT-LIGAND-PDBQT") :
                if 'Energy' in line:
                    energy.append(float(line.split("=")[1].split()[0]))
        del energy[0]
        del energy[0]
        
    return energy 

def make_pocket(protein,ligand_data,distance=5):
    protein_atom = []
    pocket = f"HEADER    {protein.split('/')[-1][:-7]}_POCKET\nCOMPND    {protein.split('/')[-1][:-7]}_POCKET\n"
    for line in open(protein).readlines():
        if line.startswith("ATOM  ") or line.startswith("HETATM"):
            record = line[0:6] # 0
            serial = int(line[6:11]) 
            atomName = line[12:16].rstrip() # 1
            altLoc = line[16] 
            resName = line[17:20].strip() #2
            chainID = line[21] #3
            if chainID == ' ':
                chainID = 'A'
            resSeq = int(line[22:26]) # 4
            x = float(line[30:38]) # 5
            y = float(line[38:46]) # 5
            z = float(line[46:54]) # 5
            occupancy = line[54:60] # 6
            tempFactor = line[60:66] # 7
            element = line[66:78]#8
            if resName != 'HOH' and not line.strip().endswith("H"):
                protein_atom.append([record, atomName,resName,chainID,resSeq,[x,y,z],occupancy,tempFactor,element])

    ligand_xyz = [data[2] for data in [*ligand_data.values()]]
    
    distance_atom = []
    inner_atom = []
    inner_chain = []
    for atom in protein_atom:
        atom_min_dis = 12
        for xyz in ligand_xyz:
            if atom_min_dis >= math.dist(atom[5],xyz):
                atom_min_dis = math.dist(atom[5],xyz)
        if atom_min_dis <= distance:
            inner_atom.append(atom[4])
            inner_chain.append(atom[3])
    inner_atom_set = set(inner_atom)
    inner_chain_set = max(set(inner_chain),key = inner_chain.count)
    serial_num = 0
    for atom in protein_atom:
        if atom[4] in inner_atom_set and atom[3] in inner_chain_set:
            serial_num += 1
            #               0:6    6:11                 12:16       17:21         21     22:26              30:38        38:46          46:54
            pocket += f"{atom[0]}{serial_num: >5} {atom[1]: <4} {atom[2]: <4}{atom[3]}{atom[4]: >4}    {atom[5][0]: >8}{atom[5][1]: >8}{atom[5][2]: >8}{atom[6]}{atom[7]}{atom[8]}\n"
    pocket += "END"
    return pocket


def pdb_list_cut(ligand_pdb):
    pdbs = []

    with open(ligand_pdb) as f:
        lines = f.read().split('ENDMDL')[:-1]
    for line in lines:
        pdbs.append(line+'ENDMDL')

    return pdbs


def make_graph(protein_path, ligand_pdb):
    lignad_name = ''
    for line in ligand_pdb.split('\n'):
        if line.startswith('COMPND'):
            ligand_name = line.split()[1]

    result = mol_2_graph(ligand_pdb,protein_path)
    node_attrs, edge_attrs, edge_indexs, error = result
    dmp = Data(x = node_attrs, edge_index = edge_indexs, edge_attr = edge_attrs, name=ligand_name)

    return dmp, error


if __name__ == "__main__" :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-r','--receptor', help='receptor pdb')
    parser.add_argument('-l','--ligands', help='ligands dlg or pdbqt or dlg, pdbqt list')
    args = parser.parse_args()
    pdb = args.receptor 
    ligands = args.ligands
    ligands = pdb_list_cut(ligands)
    for ligand in ligands:
        dmps = make_graph(pdb,ligand)
        print(dmps)

