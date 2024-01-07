import time
import pandas as pd
import numpy as np
from numpy import linalg
from numpy import dot
from numpy.linalg import norm
import os,sys
from tqdm import tqdm
import datatable as dt
from sklearn.metrics.pairwise import cosine_similarity
import pickle
import scipy.spatial 
from scipy.spatial import distance

import rdkit
from rdkit import Chem
from rdkit.Chem import GetFormalCharge, MolFromSmiles, AddHs, Draw, Descriptors, Descriptors3D
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import EmbedMolecule, MMFFOptimizeMolecule
from rdkit.Chem.rdMolDescriptors import CalcMolFormula


def similarity_search(molecule_list,filename,DT_data,data,hoips):
    d={} #dictionary to save l2 and cosine distances for list of molecules
    for mol in molecule_list:
        a={} #dictionary to save l2 and cosine distance in each iteration
        qr=hoips[mol]
        start_time = time.time()
        l2={}
        cos={}
        for i in tqdm(range(0,DT_data.shape[0])):#DT_data.shape[0]query through each molecule in compiled ~17-18million molecules DT_data.shape[0]
            q2=DT_data[i] #vector from database
            l2[data[i,0]]=np.linalg.norm(q1-q2)  #l2 calculation and asign it a dictionary
            cos[data[i,0]]=distance.cosine(q1,q2)#cosine_similarity(q1,q2[0])[0][0] #cosine calculation and asing it a dictionary 
        a['l2']=l2
        a['cosine']=cos
        d[mol]=a
    with open(filename,'wb') as hp:
        pickle.dump(d, hp)
    print('saved')
start_time = time.time()
grover_data=dt.fread('dataset/zinc_our_final_grover.csv')
DT_grover_data=grover_data[:,1:].to_numpy()
print("-dataset loading-- %5.3f seconds ---" % (time.time() - start_time))
with open('2022-Nov-22/ss_grover.pickle','rb') as fp:
    hoips_grover=pickle.load(fp)
steve=pd.read_csv('selected.csv')
steve=steve.drop_duplicates( 'smiles' )
molecule_list=steve.smiles.to_list()
filename='select_grover_suttonlab.pickle'
similarity_search(molecule_list,filename,DT_grover_data,grover_data,hoips_grover)
