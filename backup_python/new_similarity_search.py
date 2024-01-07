#!/usr/bin/env python
# coding: utf-8

# In[1]:


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


# In[2]:


import rdkit
from rdkit import Chem
from rdkit.Chem import GetFormalCharge, MolFromSmiles, AddHs, Draw, Descriptors, Descriptors3D
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import EmbedMolecule, MMFFOptimizeMolecule
from rdkit.Chem.rdMolDescriptors import CalcMolFormula


# In[4]:


start_time = time.time()
data=dt.fread('work/yatiwely/similarity_11_22_2022/dataset/zinc_our_final_rdkit_minmax.csv')
DT_data=data[:,1:].to_numpy()
print("-dataset loading-- %5.3f seconds ---" % (time.time() - start_time))


# In[5]:


with open('2022-Nov-22/ss_rdkit.pickle','rb') as fp:
    hoips=pickle.load(fp)


# In[6]:


rdkitft=['tpsa','asa','mlwt','mlp','ecc','rgy','isf','aring','ring','hdonor','hacceptor','sp3','spi','asp','hull_vol']


# In[26]:


def similarity_search(molecule_list,filename,DT_data,data,hoips):
    d={} #dictionary to save l2 and cosine distances for list of molecules
    for mol in molecule_list:
        a={} #dictionary to save l2 and cosine distance in each iteration
        qr=hoips[mol]
        q1=np.array(qr)
        start_time = time.time()
        l2={}
        cos={}
        for i in tqdm(range(0,DT_data.shape[0])):#DT_data.shape[0]   query through each molecule in compiled ~17-18million molecules DT_data.shape[0]
            q2=DT_data[i] #vector from database
            l2[data[i,0]]=np.linalg.norm(q1-q2)  #l2 calculation and asign it a dictionary
            cos[data[i,0]]=distance.cosine(q1,q2)#cosine_similarity(q1,q2[0])[0][0] #cosine calculation and asing it a dictionary 
        a['l2']=l2
        a['cosine']=cos
        d[mol]=a
    with open(filename,'wb') as hp:
        pickle.dump(d, hp)



# In[29]:


molecule_list=sorted(list(hoips.keys()))
filename='similarity_result_rdkit.pickle'
test=similarity_search(molecule_list[:10],filename,DT_data,data,aa) 


# In[54]:


start_time = time.time()
grover_data=dt.fread('/work/yatiwely/similarity_11-22-2022/dataset/zinc_our_final_grover.csv')
DT_grover_data=grover_data[:,10:].to_numpy()
print("-dataset loading-- %5.3f seconds ---" % (time.time() - start_time))
with open('2022-Nov-22/ss_grover.pickle','rb') as fp:
    hoips_grover=pickle.load(fp)


# In[55]:


molecule_list=sorted(list(hoips_grover.keys()))
filename='similarity_result_grover.pickle'


# In[ ]:


ss_grover=similarity_search(molecule_list[:1],filename,DT_grover_data,grover_data,hoips_grover) 


# In[ ]:




