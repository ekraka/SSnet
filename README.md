# SSnet
SSnet - Secondary Structure based End-to-End Learning model for Protein-Ligand Interaction Prediction


## Training

ECF is the morgan fingerprint for ligands (more details in the paper).

GNN is graph based neural network for ligands.

Click into corresponding directory and read the description in README.md file for further processing


## GRAD-CAM

The grad-cam is computed through the weights optimized from DUD-E dataset

### Preprocessing

Change the weight_file in the last few lines in GRAD/grad_cam.py to point GRAD/multiBranch_Conv_model_copy.h5. Check GRAD/grad_cam.py

### Usage

python path_to_GRAD/grad_cam.py pdb_file.pdb 'smiles'

Please make sure to capitalize the four letter pdb id

### Example

python /users/nirajv/PDiML/GRAD/grad_cam.py 6M18.pdb 'C[C@H](N[C@@H](CCc1ccccc1)C(=O)O)C(=O)N1[C@H](C(=O)O)C[C@H]2CCCC[C@@H]21'



