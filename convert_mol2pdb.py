#!/usr/bin/env python

from openbabel.openbabel import *

def mol2pdb(file,output):
    if not output.endswith('.pdb'):
        output += '.pdb'
    fs = open(output,'w')
    obconversion = OBConversion()
    obconversion.SetInAndOutFormats("mol2","pdb")
    obmol = OBMol()

    notatend = obconversion.ReadFile(obmol,file)
    while notatend:
        fs.write(obconversion.WriteString(obmol))
        fs.write('ENDMDL\n')
        obmol = OBMol()
        notatend = obconversion.Read(obmol)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-l','--ligands', default='examples/1nc1_ligand.mol2', help='ligands .mol2')
    parser.add_argument('-o','--output', default='examples/1nc1_mol2_ligand.pdb', help='output .pdb')
    args = parser.parse_args()
    
    mol2pdb(args.ligands,args.output)
