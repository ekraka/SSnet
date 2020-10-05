
### Usage

parameters:
-m <mode>
    > mode = 10 (Default)
        (This model was trained based on 10nM IC50 cutoff for actives)
    > mode = 25
        (This model was trained based on 25nM IC50 cutoff for actives)
    > mode = 100
        (This model was trained based on 100nM IC50 cutoff for actives)
    > mode = grad
        (This model was trained for generating heatmaps for potential binding location)
  
 -h
    provides help I guess!


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
The input file should math the format as given in 'example_multi_input.txt'
i.e smile [space] target


### Outputs

results_pdb_file.txt

### Requirements (python modules)

1) wget
2) keras
3) tensorflow
4) rdkit





## GRAD CAM analysis

GradCAM analysis can be perfomed by changing mode to grad
-m grad

### Example

python path_to_SSnet_in_Action/predict.py -t 6M18.pdb -l 'C[C@H](N[C@@H](CCc1ccccc1)C(=O)O)C(=O)N1[C@H](C(=O)O)C[C@H]2CCCC[C@@H]21' -m grad

### Outputs

results_pdb_file.pdb

The heatmap is wriiten as b-factor

To visualize in pymol type the following command after loading the pdb:

spectrum b, rainbow, minimum = 0.0, maximum = 100.0


