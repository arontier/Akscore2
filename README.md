# Akscore2


# input data
```
#단일
python input_data.py -r input_example/1a30_protein.pdb -l input_example/test_pdb.pdb
#멀티
python input_data.py -r input_example/1a30_protein.pdb -l input_example/1a30_ligand.in.pdb

pdb_list_cut()으로 ligand pdb파일을 받아주고, output을 for문 돌려서 make_graph에 protein_pdb와 같이 받아주면 된다.

```

# input error

```
python input_data.py -r input_example/1a30_protein.pdb -l input_example/test_pdb_error.pdb
```
첫번째 ligand만 에러남. error log로 error ligand가 뜸
