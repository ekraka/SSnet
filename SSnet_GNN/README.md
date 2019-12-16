

## Usage

We provide two major scripts:

- code/preprocess_data.py creates the input tensor data of CPIs
for processing with PyTorch from the original data
(see dataset/human or celegans/original/data.txt).
- code/run_training.py trains the model using the above preprocessed data
(see dataset/human or celegans/input).

(i) Create the tensor data of CPIs with the following command:
```
cd code
bash preprocess_data.sh
```

The preprocessed data are saved in the dataset/input directory.

(ii) Using the preprocessed data, train the model with the following command:
```
bash run_training.sh
```

The training and test results and the model are saved in the output directory
(after training, see output/result and output/model).

(iii) You can change the hyperparameters in preprocess_data.sh and run_training.sh.
Try to learn various models.


## Result

Learning curves (x-axis is epoch and y-axis is AUC)
on the test datasets of human and *C. elegans* are as follows:

<div align="center">
<p><img src="learning_curves.jpeg" width="800" /></p>
</div>

These results can be reproduce by the above two commands (i) and (ii).
