#!/bin/sh
#SBATCH --job-name=001
#SBATCH --output 001.out
#SBATCH --error 001.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p suttonlab-48core

module load python3/anaconda/2021.07
source activate /work/yatiwely/ENVS/rdkit
hostname
date

python3 dictionary_creation.py 001.csv
