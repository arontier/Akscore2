# Akscore2


### input data
```
#단일
python input_data.py -r input_example/1a30_protein.pdb -l input_example/test_pdb.pdb
#멀티
python input_data.py -r input_example/1a30_protein.pdb -l input_example/1a30_ligand.in.pdb

pdb_list_cut()으로 ligand pdb파일을 받아주고, output을 for문 돌려서 make_graph에 protein_pdb와 같이 받아주면 된다.

```

multi output
```
(Data(x=[127, 73], edge_index=[2, 2100], edge_attr=[2100, 24], name='1a30_128'), 1)
(Data(x=[165, 73], edge_index=[2, 2126], edge_attr=[2126, 24], name='1a30_130'), 1)
(Data(x=[164, 73], edge_index=[2, 2150], edge_attr=[2150, 24], name='1a30_139'), 1)
(Data(x=[117, 73], edge_index=[2, 1720], edge_attr=[1720, 24], name='1a30_140'), 1)
(Data(x=[172, 73], edge_index=[2, 2236], edge_attr=[2236, 24], name='1a30_144'), 1)
(Data(x=[142, 73], edge_index=[2, 1696], edge_attr=[1696, 24], name='1a30_167'), 1)
(Data(x=[139, 73], edge_index=[2, 1940], edge_attr=[1940, 24], name='1a30_170'), 1)
(Data(x=[147, 73], edge_index=[2, 1890], edge_attr=[1890, 24], name='1a30_186'), 1)
(Data(x=[133, 73], edge_index=[2, 1704], edge_attr=[1704, 24], name='1a30_187'), 1)
(Data(x=[153, 73], edge_index=[2, 2192], edge_attr=[2192, 24], name='1a30_221'), 1)
(Data(x=[133, 73], edge_index=[2, 1716], edge_attr=[1716, 24], name='1a30_222'), 1)
(Data(x=[157, 73], edge_index=[2, 1642], edge_attr=[1642, 24], name='1a30_22'), 1)
(Data(x=[154, 73], edge_index=[2, 2076], edge_attr=[2076, 24], name='1a30_267'), 1)
(Data(x=[150, 73], edge_index=[2, 2004], edge_attr=[2004, 24], name='1a30_270'), 1)
(Data(x=[163, 73], edge_index=[2, 1508], edge_attr=[1508, 24], name='1a30_284'), 1)
(Data(x=[148, 73], edge_index=[2, 1862], edge_attr=[1862, 24], name='1a30_296'), 1)
(Data(x=[173, 73], edge_index=[2, 2456], edge_attr=[2456, 24], name='1a30_314'), 1)
(Data(x=[138, 73], edge_index=[2, 2152], edge_attr=[2152, 24], name='1a30_319'), 1)
(Data(x=[145, 73], edge_index=[2, 2084], edge_attr=[2084, 24], name='1a30_336'), 1)
(Data(x=[146, 73], edge_index=[2, 2016], edge_attr=[2016, 24], name='1a30_348'), 1)
(Data(x=[138, 73], edge_index=[2, 1990], edge_attr=[1990, 24], name='1a30_353'), 1)
(Data(x=[131, 73], edge_index=[2, 2142], edge_attr=[2142, 24], name='1a30_35'), 1)
(Data(x=[159, 73], edge_index=[2, 2148], edge_attr=[2148, 24], name='1a30_388'), 1)
(Data(x=[135, 73], edge_index=[2, 2062], edge_attr=[2062, 24], name='1a30_407'), 1)
(Data(x=[119, 73], edge_index=[2, 2024], edge_attr=[2024, 24], name='1a30_418'), 1)
(Data(x=[127, 73], edge_index=[2, 2190], edge_attr=[2190, 24], name='1a30_450'), 1)
(Data(x=[124, 73], edge_index=[2, 1364], edge_attr=[1364, 24], name='1a30_452'), 1)
(Data(x=[145, 73], edge_index=[2, 2158], edge_attr=[2158, 24], name='1a30_490'), 1)
(Data(x=[127, 73], edge_index=[2, 2258], edge_attr=[2258, 24], name='1a30_498'), 1)
(Data(x=[120, 73], edge_index=[2, 2018], edge_attr=[2018, 24], name='1a30_538'), 1)
(Data(x=[141, 73], edge_index=[2, 1904], edge_attr=[1904, 24], name='1a30_569'), 1)
(Data(x=[139, 73], edge_index=[2, 1528], edge_attr=[1528, 24], name='1a30_664'), 1)
(Data(x=[112, 73], edge_index=[2, 1634], edge_attr=[1634, 24], name='1a30_670'), 1)
(Data(x=[154, 73], edge_index=[2, 1646], edge_attr=[1646, 24], name='1a30_671'), 1)
(Data(x=[131, 73], edge_index=[2, 2102], edge_attr=[2102, 24], name='1a30_6'), 1)
(Data(x=[157, 73], edge_index=[2, 2070], edge_attr=[2070, 24], name='1a30_704'), 1)
(Data(x=[153, 73], edge_index=[2, 2040], edge_attr=[2040, 24], name='1a30_705'), 1)
(Data(x=[153, 73], edge_index=[2, 2036], edge_attr=[2036, 24], name='1a30_706'), 1)
(Data(x=[157, 73], edge_index=[2, 2072], edge_attr=[2072, 24], name='1a30_709'), 1)
(Data(x=[157, 73], edge_index=[2, 2072], edge_attr=[2072, 24], name='1a30_712'), 1)
(Data(x=[157, 73], edge_index=[2, 2086], edge_attr=[2086, 24], name='1a30_718'), 1)
(Data(x=[157, 73], edge_index=[2, 1986], edge_attr=[1986, 24], name='1a30_731'), 1)
(Data(x=[134, 73], edge_index=[2, 1608], edge_attr=[1608, 24], name='1a30_733'), 1)
(Data(x=[131, 73], edge_index=[2, 1726], edge_attr=[1726, 24], name='1a30_736'), 1)
(Data(x=[157, 73], edge_index=[2, 2106], edge_attr=[2106, 24], name='1a30_743'), 1)
(Data(x=[123, 73], edge_index=[2, 1610], edge_attr=[1610, 24], name='1a30_744'), 1)
(Data(x=[146, 73], edge_index=[2, 2062], edge_attr=[2062, 24], name='1a30_747'), 1)
(Data(x=[127, 73], edge_index=[2, 1636], edge_attr=[1636, 24], name='1a30_754'), 1)
(Data(x=[133, 73], edge_index=[2, 1828], edge_attr=[1828, 24], name='1a30_759'), 1)
(Data(x=[142, 73], edge_index=[2, 2034], edge_attr=[2034, 24], name='1a30_75'), 1)
(Data(x=[146, 73], edge_index=[2, 2078], edge_attr=[2078, 24], name='1a30_797'), 1)
(Data(x=[149, 73], edge_index=[2, 2024], edge_attr=[2024, 24], name='1a30_828'), 1)
(Data(x=[135, 73], edge_index=[2, 1776], edge_attr=[1776, 24], name='1a30_868'), 1)
(Data(x=[138, 73], edge_index=[2, 2218], edge_attr=[2218, 24], name='1a30_86'), 1)
(Data(x=[117, 73], edge_index=[2, 2094], edge_attr=[2094, 24], name='1a30_87'), 1)
(Data(x=[120, 73], edge_index=[2, 1892], edge_attr=[1892, 24], name='1a30_883'), 1)
(Data(x=[130, 73], edge_index=[2, 1948], edge_attr=[1948, 24], name='1a30_899'), 1)
(Data(x=[134, 73], edge_index=[2, 2256], edge_attr=[2256, 24], name='1a30_90'), 1)
(Data(x=[142, 73], edge_index=[2, 2056], edge_attr=[2056, 24], name='1a30_920'), 1)
(Data(x=[143, 73], edge_index=[2, 2076], edge_attr=[2076, 24], name='1a30_935'), 1)
(Data(x=[143, 73], edge_index=[2, 2008], edge_attr=[2008, 24], name='1a30_938'), 1)
(Data(x=[134, 73], edge_index=[2, 1666], edge_attr=[1666, 24], name='1a30_940'), 1)
(Data(x=[142, 73], edge_index=[2, 1642], edge_attr=[1642, 24], name='1a30_942'), 1)
(Data(x=[135, 73], edge_index=[2, 1982], edge_attr=[1982, 24], name='1a30_944'), 1)
(Data(x=[135, 73], edge_index=[2, 2074], edge_attr=[2074, 24], name='1a30_955'), 1)
(Data(x=[135, 73], edge_index=[2, 1956], edge_attr=[1956, 24], name='1a30_972'), 1)
(Data(x=[135, 73], edge_index=[2, 2086], edge_attr=[2086, 24], name='1a30_985'), 1)
(Data(x=[149, 73], edge_index=[2, 2178], edge_attr=[2178, 24], name='1a30_992'), 1)
(Data(x=[131, 73], edge_index=[2, 2100], edge_attr=[2100, 24], name='1a30_997'), 1)
(Data(x=[138, 73], edge_index=[2, 2086], edge_attr=[2086, 24], name='1a30_999'), 1)
(Data(x=[153, 73], edge_index=[2, 2014], edge_attr=[2014, 24], name='1a30_ligand'), 1)
```
### input error

```
python input_data.py -r input_example/1a30_protein.pdb -l input_example/test_pdb_error.pdb
```
첫번째 ligand만 에러남. error log로 error ligand가 뜸

output
```
INFO:root:, error ligand
(Data(x=[100, 73], edge_index=[2, 1000], edge_attr=[1000, 24], name='1a30_128'), 0)
(Data(x=[165, 73], edge_index=[2, 2126], edge_attr=[2126, 24], name='1a30_130'), 1)
(Data(x=[164, 73], edge_index=[2, 2150], edge_attr=[2150, 24], name='1a30_139'), 1)
(Data(x=[117, 73], edge_index=[2, 1720], edge_attr=[1720, 24], name='1a30_140'), 1)
(Data(x=[172, 73], edge_index=[2, 2236], edge_attr=[2236, 24], name='1a30_144'), 1)
(Data(x=[142, 73], edge_index=[2, 1696], edge_attr=[1696, 24], name='1a30_167'), 1)
(Data(x=[139, 73], edge_index=[2, 1940], edge_attr=[1940, 24], name='1a30_170'), 1)
(Data(x=[147, 73], edge_index=[2, 1890], edge_attr=[1890, 24], name='1a30_186'), 1)
(Data(x=[133, 73], edge_index=[2, 1704], edge_attr=[1704, 24], name='1a30_187'), 1)
(Data(x=[153, 73], edge_index=[2, 2192], edge_attr=[2192, 24], name='1a30_221'), 1)
(Data(x=[133, 73], edge_index=[2, 1716], edge_attr=[1716, 24], name='1a30_222'), 1)
(Data(x=[157, 73], edge_index=[2, 1642], edge_attr=[1642, 24], name='1a30_22'), 1)
(Data(x=[154, 73], edge_index=[2, 2076], edge_attr=[2076, 24], name='1a30_267'), 1)
(Data(x=[150, 73], edge_index=[2, 2004], edge_attr=[2004, 24], name='1a30_270'), 1)
(Data(x=[163, 73], edge_index=[2, 1508], edge_attr=[1508, 24], name='1a30_284'), 1)
(Data(x=[148, 73], edge_index=[2, 1862], edge_attr=[1862, 24], name='1a30_296'), 1)
(Data(x=[173, 73], edge_index=[2, 2456], edge_attr=[2456, 24], name='1a30_314'), 1)
(Data(x=[138, 73], edge_index=[2, 2152], edge_attr=[2152, 24], name='1a30_319'), 1)
(Data(x=[145, 73], edge_index=[2, 2084], edge_attr=[2084, 24], name='1a30_336'), 1)
(Data(x=[146, 73], edge_index=[2, 2016], edge_attr=[2016, 24], name='1a30_348'), 1)
(Data(x=[138, 73], edge_index=[2, 1990], edge_attr=[1990, 24], name='1a30_353'), 1)
(Data(x=[131, 73], edge_index=[2, 2142], edge_attr=[2142, 24], name='1a30_35'), 1)
(Data(x=[159, 73], edge_index=[2, 2148], edge_attr=[2148, 24], name='1a30_388'), 1)
(Data(x=[135, 73], edge_index=[2, 2062], edge_attr=[2062, 24], name='1a30_407'), 1)
(Data(x=[119, 73], edge_index=[2, 2024], edge_attr=[2024, 24], name='1a30_418'), 1)
(Data(x=[127, 73], edge_index=[2, 2190], edge_attr=[2190, 24], name='1a30_450'), 1)
(Data(x=[124, 73], edge_index=[2, 1364], edge_attr=[1364, 24], name='1a30_452'), 1)
(Data(x=[145, 73], edge_index=[2, 2158], edge_attr=[2158, 24], name='1a30_490'), 1)
(Data(x=[127, 73], edge_index=[2, 2258], edge_attr=[2258, 24], name='1a30_498'), 1)
(Data(x=[120, 73], edge_index=[2, 2018], edge_attr=[2018, 24], name='1a30_538'), 1)
(Data(x=[141, 73], edge_index=[2, 1904], edge_attr=[1904, 24], name='1a30_569'), 1)
(Data(x=[139, 73], edge_index=[2, 1528], edge_attr=[1528, 24], name='1a30_664'), 1)
(Data(x=[112, 73], edge_index=[2, 1634], edge_attr=[1634, 24], name='1a30_670'), 1)
(Data(x=[154, 73], edge_index=[2, 1646], edge_attr=[1646, 24], name='1a30_671'), 1)
(Data(x=[131, 73], edge_index=[2, 2102], edge_attr=[2102, 24], name='1a30_6'), 1)
(Data(x=[157, 73], edge_index=[2, 2070], edge_attr=[2070, 24], name='1a30_704'), 1)
(Data(x=[153, 73], edge_index=[2, 2040], edge_attr=[2040, 24], name='1a30_705'), 1)
(Data(x=[153, 73], edge_index=[2, 2036], edge_attr=[2036, 24], name='1a30_706'), 1)
(Data(x=[157, 73], edge_index=[2, 2072], edge_attr=[2072, 24], name='1a30_709'), 1)
(Data(x=[157, 73], edge_index=[2, 2072], edge_attr=[2072, 24], name='1a30_712'), 1)
(Data(x=[157, 73], edge_index=[2, 2086], edge_attr=[2086, 24], name='1a30_718'), 1)
(Data(x=[157, 73], edge_index=[2, 1986], edge_attr=[1986, 24], name='1a30_731'), 1)
(Data(x=[134, 73], edge_index=[2, 1608], edge_attr=[1608, 24], name='1a30_733'), 1)
(Data(x=[131, 73], edge_index=[2, 1726], edge_attr=[1726, 24], name='1a30_736'), 1)
(Data(x=[157, 73], edge_index=[2, 2106], edge_attr=[2106, 24], name='1a30_743'), 1)
(Data(x=[123, 73], edge_index=[2, 1610], edge_attr=[1610, 24], name='1a30_744'), 1)
(Data(x=[146, 73], edge_index=[2, 2062], edge_attr=[2062, 24], name='1a30_747'), 1)
(Data(x=[127, 73], edge_index=[2, 1636], edge_attr=[1636, 24], name='1a30_754'), 1)
(Data(x=[133, 73], edge_index=[2, 1828], edge_attr=[1828, 24], name='1a30_759'), 1)
(Data(x=[142, 73], edge_index=[2, 2034], edge_attr=[2034, 24], name='1a30_75'), 1)
(Data(x=[146, 73], edge_index=[2, 2078], edge_attr=[2078, 24], name='1a30_797'), 1)
(Data(x=[149, 73], edge_index=[2, 2024], edge_attr=[2024, 24], name='1a30_828'), 1)
(Data(x=[135, 73], edge_index=[2, 1776], edge_attr=[1776, 24], name='1a30_868'), 1)
(Data(x=[138, 73], edge_index=[2, 2218], edge_attr=[2218, 24], name='1a30_86'), 1)
(Data(x=[117, 73], edge_index=[2, 2094], edge_attr=[2094, 24], name='1a30_87'), 1)
(Data(x=[120, 73], edge_index=[2, 1892], edge_attr=[1892, 24], name='1a30_883'), 1)
(Data(x=[130, 73], edge_index=[2, 1948], edge_attr=[1948, 24], name='1a30_899'), 1)
(Data(x=[134, 73], edge_index=[2, 2256], edge_attr=[2256, 24], name='1a30_90'), 1)
(Data(x=[142, 73], edge_index=[2, 2056], edge_attr=[2056, 24], name='1a30_920'), 1)
(Data(x=[143, 73], edge_index=[2, 2076], edge_attr=[2076, 24], name='1a30_935'), 1)
(Data(x=[143, 73], edge_index=[2, 2008], edge_attr=[2008, 24], name='1a30_938'), 1)
(Data(x=[134, 73], edge_index=[2, 1666], edge_attr=[1666, 24], name='1a30_940'), 1)
(Data(x=[142, 73], edge_index=[2, 1642], edge_attr=[1642, 24], name='1a30_942'), 1)
(Data(x=[135, 73], edge_index=[2, 1982], edge_attr=[1982, 24], name='1a30_944'), 1)
(Data(x=[135, 73], edge_index=[2, 2074], edge_attr=[2074, 24], name='1a30_955'), 1)
(Data(x=[135, 73], edge_index=[2, 1956], edge_attr=[1956, 24], name='1a30_972'), 1)
(Data(x=[135, 73], edge_index=[2, 2086], edge_attr=[2086, 24], name='1a30_985'), 1)
(Data(x=[149, 73], edge_index=[2, 2178], edge_attr=[2178, 24], name='1a30_992'), 1)
(Data(x=[131, 73], edge_index=[2, 2100], edge_attr=[2100, 24], name='1a30_997'), 1)
(Data(x=[138, 73], edge_index=[2, 2086], edge_attr=[2086, 24], name='1a30_999'), 1)
(Data(x=[153, 73], edge_index=[2, 2014], edge_attr=[2014, 24], name='1a30_ligand'), 1)

```
## License
Shield: [![CC BY-NC-SA 4.0][cc-by-nc-sa-shield]][cc-by-nc-sa]

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/
[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png
[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg
