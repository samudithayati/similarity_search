def bash_file_creater(filename):
    xyz =open(filename.split('.')[0]+'.sh','w')
    xyz.writelines('#!/bin/sh\n')
    xyz.writelines('#SBATCH --job-name='+filename.split('.')[0]+'\n')
    xyz.writelines('#SBATCH --output '+filename.split('.')[0]+'.out\n')
    xyz.writelines('#SBATCH --error '+filename.split('.')[0]+'.err\n')
    xyz.writelines('#SBATCH -N 1\n')
    xyz.writelines('#SBATCH -n 1\n')
    xyz.writelines('#SBATCH -p suttonlab-48core\n')
    xyz.writelines('\n')
    xyz.writelines('module load python3/anaconda/2021.07\n')
    xyz.writelines('source activate /work/yatiwely/ENVS/rdkit\n')
    xyz.writelines('hostname\n')
    xyz.writelines('date\n')
    xyz.writelines('\n')
    
    xyz.writelines('python3 dictionary_creation.py '+filename+'\n')
    xyz.close()    
    
import glob
files=glob.glob('*.csv')
for ii in files:
    bash_file_creater(ii)
