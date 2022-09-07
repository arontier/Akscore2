# Akscore2
## Requirements
- python 
- pytorch == 1.11.0 (cuda==11.3, cudnn==8.2.0)
- pytorch geometric == 2.0.4
- rdkit == 2022.3.1
- openbabel == 3.1.1
- meeko == 0.3.3
- numpy >= 1.22
- pandas >= 1.4.2

## Usages
(1) Convert ligand .mol2 files into .pdb format that like in the example directory.
```
python convert_mol2pdb.py -l examples/1nc1_ligand.mol2 -o examples/1nc1_mol2_ligand.pdb
```
(2) Run Akscore2 to get binding affinity score by following command. It will first construct graph using .pdb file and then output result of Akscore2_NonDock, Akscore2_Dock and their multiplication. 
```
python akscore2_run.py -r examples/1nc1_protein.pdb -l examples/1nc1_ligand.pdb -s akscore2_dockc -o result.csv --batch_size 8 --num_workers 8 --device gpu
```
If you want to use Akscore2_DockS model instead of Akscore2_DockC, change to `-s akscore2_docks`.

In the result file, those failed to contruct graph are filled with NaN.

(3) If you want to get final ensemble score described in the paper that includes autodock score, you need to install and run the autodock by yourself and simply plus the score to the result of Akscore2. 

## License
Shield: [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg
