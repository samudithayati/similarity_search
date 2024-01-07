#!/bin/sh
#SBATCH --job-name=l2-chiral2
#SBATCH --output bash_out/l2_grover%j.out
#SBATCH --error bash_out/l2_grover%j.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p nodename
##SBATCH --mail-type=FAIL,END
##SBATCH --mail-user=

##Load your modules first:

module load python3/anaconda/2021.07
source activate /work/username/ENVS/rdkit
##Add your code here:

hostname
date
                 
python run_similarity.py
