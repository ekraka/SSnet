
### Usage

python path_to_SSnet_in_Action/predict.py pdb_file.pdb 'smiles'

Please make sure to capitalize the four letter pdb id if the pdb is not present in the working directory (SSnet will download from RCSB website if present).

### Example

python path_to_SSnet_in_Action/predict.py 6M18.pdb 'C[C@H](N[C@@H](CCc1ccccc1)C(=O)O)C(=O)N1[C@H](C(=O)O)C[C@H]2CCCC[C@@H]21'


#### Provide ligands as .smi file (eg test.smi)

python path_to_SSnet_in_Action/predict.py 6M18.pdb test.smi

### Outputs

results_pdb_file.txt

### Requirements (python modules)

1) wget
2) keras
3) tensorflow
4) rdkit

