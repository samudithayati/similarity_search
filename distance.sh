#!/bin/sh
#SBATCH --job-name=l2-chiral2
#SBATCH --output bash_out/l2_grover%j.out
#SBATCH --error bash_out/l2_grover%j.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p BigMem
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=yatiwely@email.sc.edu

##Load your modules first:

module load python3/anaconda/2021.07
source activate /work/yatiwely/ENVS/rdkit
##Add your code here:

hostname
date
                 
python script/grover1.py
