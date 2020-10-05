
### Usage
## Screening
Screening can be performed based on single target - single ligand, single target - multiple ligand as: 

A protein is provided by the parameter '-t',
Ligand(s) are provided as '-l' 

python path_to_SSnet_in_Action/predict.py -t pdb_file.pdb -l 'smiles'

Please make sure to capitalize the four letter pdb id if the protein is not present in the working directory (SSnet will download from RCSB website if present).

### Example

python path_to_SSnet_in_Action/predict.py -t 6M18.pdb -l 'C[C@H](N[C@@H](CCc1ccccc1)C(=O)O)C(=O)N1[C@H](C(=O)O)C[C@H]2CCCC[C@@H]21'


#### Provide ligands as .smi file (eg test.smi) for multiple ligands

python path_to_SSnet_in_Action/predict.py -t 6M18.pdb -l test.smi

#### For multiple target and multiple ligands

use '-i' followed by input file


### Outputs

results_pdb_file.txt

### Requirements (python modules)

1) wget
2) keras
3) tensorflow
4) rdkit

