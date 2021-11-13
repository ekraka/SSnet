
### Usage

### Binary files for Linux, Windows and OS-X can be found at https://smu.box.com/s/yjd4dm1k9c4axmg5i0wihbz3nx3tx1er
Note that binary files do not require any pre-installation such as python. Replace 'python path_to_SSnet_in_Action/predict.py' to the binary file. Check below for arguments.

If using python, make sure you have the following modules:

    keras
    h5py 2.10
    python-wget
    rdkit

parameters:
Note: Based on our experience, -m 25 works best for enrichment!


-m <mode>
    
      mode = 10 (Default)
        (This model was trained based on 10nM IC50 cutoff for actives)
    
      mode = 25
        (This model was trained based on 25nM IC50 cutoff for actives)
        
      mode = 100
        (This model was trained based on 100nM IC50 cutoff for actives)
        
      mode = grad
        (This model was trained for generating heatmaps for potential binding location)
  
 -h
    provides help I guess!
    
        python path_to_SSnet_in_Action/predict.py -h


## Screening
Screening can be performed based on single target - single ligand, single target - multiple ligand as: 

A protein is provided by the parameter '-t',
Ligand(s) are provided as '-l' 

    python path_to_SSnet_in_Action/predict.py -t pdb_file.pdb -l 'smiles'

Please make sure to capitalize the four letter pdb id if the protein is not present in the working directory (SSnet will download from RCSB website if present).

### Example

    python path_to_SSnet_in_Action/predict.py -t 6M18.pdb -l 'C[C@H](N[C@@H](CCc1ccccc1)C(=O)O)C(=O)N1[C@H](C(=O)O)C[C@H]2CCCC[C@@H]21'
    
    python path_to_SSnet_in_Action/predict.py -m 25 -t 6M18.pdb -l 'C[C@H](N[C@@H](CCc1ccccc1)C(=O)O)C(=O)N1[C@H](C(=O)O)C[C@H]2CCCC[C@@H]21'


#### Provide ligands as .smi file (eg test.smi) for multiple ligands

    python path_to_SSnet_in_Action/predict.py -t 6M18.pdb -l test.smi

#### For multiple target and multiple ligands

use '-i' followed by input file.

The input file should match the format as given in 'example_multi_input.txt'

i.e smile [space] target

-t_dir points to the directory where pdb files are (if not present, SSnet will try to download from RCSB website)

    python path_to_SSnet_in_Action/predict.py -i example_multi_input.txt -t_dir '.'
    
In the above case, -t_dir points to the current working directory


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
    
Note that GradCAM is ligand specific and thus only takes SMILES as ligand input (Multiple smiles through .smi files cannot be executed)

### Example

    python path_to_SSnet_in_Action/predict.py -t 6M18.pdb -l 'C[C@H](N[C@@H](CCc1ccccc1)C(=O)O)C(=O)N1[C@H](C(=O)O)C[C@H]2CCCC[C@@H]21' -m grad

### Outputs

pdb_file_GRAD.pdb

The heatmap is written as b-factor

To visualize in pymol type the following command after loading the pdb:

    spectrum b, rainbow, minimum = 0.0, maximum = 100.0


