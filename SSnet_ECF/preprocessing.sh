#!/bin/bash

inps=inp 

python preprocessing_ligand.py $inps

python preprocessing_target.py $inps

python rdk.py $inps

