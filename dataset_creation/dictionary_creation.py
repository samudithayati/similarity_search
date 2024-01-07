#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd
import glob
import sys,os,pickle
#convex hull pacakge import 
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#rdkit import 
import rdkit
from rdkit import Chem
from rdkit.Chem import GetFormalCharge, MolFromSmiles, AddHs, Draw, Descriptors, Descriptors3D,rdMolDescriptors
from rdkit.Chem import AllChem , Draw
from rdkit.Chem.AllChem import EmbedMolecule, MMFFOptimizeMolecule
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

import json
from tqdm import tqdm
import seaborn as sns
import datetime
import random


# In[2]:


random.seed(42)


# In[3]:


def delete_file(filename):
    if os.path.exists(filename):
        os.remove(filename)
        #print("The file has been deleted successfully")
    else:
        pass
        #print("The file does not exist!")

def create_directory(dir_name):
    if os.path.exists(dir_name):
        pass
        #print('previous path exist')
    else:
        os.makedirs(dir_name)
    return(dir_name)


# In[4]:


def mol_filereader(fn):
    x,y,z,atom = ([] for i in range(4))
    xyz = open(fn, "r")
    xyz.readline()
    xyz.readline()
    xyz.readline()
    to=xyz.readline().split()
    total_atom=int(to[0])
    for i in range(0,total_atom):
        line_data=xyz.readline().split()
        x.append(float(line_data[0]))
        y.append(float(line_data[1]))
        z.append(float(line_data[2]))
        atom.append(line_data[3])
    return(x,y,z,atom)

def smiles_xyz(m3,fname):

    file = open(x_file+"/"+str(fname)+'.mol','w+')                
    file.write(Chem.MolToMolBlock(m3)) 
    file.close()
    x,y,z,atom=mol_filereader(x_file+"/"+str(fname)+'.mol')
    xyz =open(x_file+"/"+str(fname)+'.xyz','w')
    xyz.writelines(str(len(x))+'\n')
    xyz.writelines('\n')
    for ii in range(0,len(x)):
        xyz.writelines(str(atom[ii])+' '+str(round(x[ii],2))+' '+str(round(y[ii],2))+' '+str(round(z[ii],2))+'\n')
    xyz.writelines('\n')
    xyz.close()
    return(x,y,z)

def hull_volume(x,y,z):
    coord=[]
    for i in range(0,len(x)):
        cd=[float(x[i]),float(y[i]),float(z[i])]
        coord.append(cd)
    res= np.array(coord)
    try:
        hull = ConvexHull(res)
        return(hull.volume)
    except:
        return(0)


# In[5]:


def convert_smiles_to_canonical(smile):
#     print(smile)
    try:
        up_smile=Chem.MolToSmiles(Chem.MolFromSmiles(smile))
    except:
        up_smile=np.nan
    return(up_smile)


# In[6]:


filename=sys.argv[1]
x_file=filename.split('.')[0]
create_directory(x_file)


# In[7]:



unfiltered_data=pd.read_csv(filename)


# In[8]:



tmp_smiles=unfiltered_data.smiles.to_list()
print(len(tmp_smiles))
tmp_smiles=[x for x in tmp_smiles if str(x) != 'nan'] #remove nan
print(len(tmp_smiles))


unique_cations=sorted(set(tmp_smiles))
print(len(unique_cations))


# In[9]:


def feature_gen_xyz(SMILES):
# SMILES=unique_cations[:5]
    def tpsa_f(ml):return(Descriptors.TPSA(ml))
    def asa_f(ml):return(Descriptors.LabuteASA(ml))
    # def hd_f(ml):return(Descriptors.NumHDonors(ml)) #same as rdMolDescriptors.CalcNumHBD
    def mlwt_f(ml):return(Descriptors.MolWt(ml))
    def mlp_f(ml):return(Descriptors.MolLogP(ml))
    def ecc_f(ml):return(Descriptors3D.Eccentricity(ml))
    def rgy_f(ml):return(Descriptors3D.RadiusOfGyration(ml))
    def isf_f(ml):return(Descriptors3D.InertialShapeFactor(ml))
    def pm1_f(ml):return(Descriptors3D.PMI1(ml))
    def pm2_f(ml):return(Descriptors3D.PMI2(ml))
    def pm3_f(ml):return(Descriptors3D.PMI3(ml))
    def aring_f(ml):return(rdMolDescriptors.CalcNumAromaticRings(ml))
    def ring_f(ml):return(rdMolDescriptors.CalcNumRings(ml))
    def hdonor_f(ml):return(rdMolDescriptors.CalcNumHBD(ml))
    def hacceptor_f(ml):return(rdMolDescriptors.CalcNumHBA(ml))
    def sp3_f(ml):return(rdMolDescriptors.CalcFractionCSP3(ml))
    def spi_f(ml):return(Descriptors3D.SpherocityIndex(ml))
    def asp_f(ml):return(Descriptors3D.Asphericity(ml))
    def hull_f(ml):return()
    
    list_fun=[tpsa_f,asa_f,mlwt_f,mlp_f,ecc_f,rgy_f,isf_f,pm1_f,pm2_f,pm3_f,aring_f,ring_f,hdonor_f,hacceptor_f,sp3_f,spi_f,asp_f]
    list_prop=['tpsa','asa','mlwt','mlp','ecc','rgy','isf','pm1','pm2','pm3','aring','ring','hdonor','hacceptor','sp3','spi','asp']

    function_dict={}
    for i1,i2 in zip(list_prop,list_fun):
        function_dict[i1]=i2

    def function(abc,p,ml):
        try:
            prop_val=abc[p](ml)
        except:
            prop_val=np.nan
        return(prop_val)

    smile=[]
    function_dict={}
    for i1,i2 in zip(list_prop,list_fun):
        function_dict[i1]=i2

    invalid=[]

    molecule_dict={}
    xyz_name={}
    random_seed=42
    for ii in tqdm(enumerate(SMILES)):
        try:
            mol=MolFromSmiles(ii[1])
            molh=AddHs(mol)
            smile.append(ii[1])
            try:
                Chem.AllChem.EmbedMolecule(molh,maxAttempts=500,randomSeed=42)#,useRandomCoords=True,)
                Chem.rdForceFieldHelpers.MMFFOptimizeMolecule(molh)
                AllChem.MMFFOptimizeMolecule(molh, maxIters=500)
                tmp_dict={}
      
                for i in list_prop:
                    tmp_dict[i]=function(function_dict,i,molh)

                x,y,z=smiles_xyz(molh,ii[0])
#                 xyz_name[ii[1]]=ii[0]
#                 print(hull_volume(x,y,z))
                try:
                    tmp_dict['hull_vol']=hull_volume(x,y,z)
                    
                except:
                    tmp_dict['hull_vol']=np.nan
                tmp_dict['filename']=filename.split('.')[0]+str(ii[0])
#                 molecule_dict[ii[1]]=tmp_dict
            except:
                tmp_dict={}
                for i in list_prop:
                    tmp_dict[i]=np.nan
                tmp_dict['hull_vol']=np.nan
                tmp_dict['filename']=filename.split('.')[0]+str(ii[0])
            molecule_dict[ii[1]]=tmp_dict

        except:
            invalid.append(ii[1])

    return(molecule_dict)


# In[10]:


unique_cation_dict= feature_gen_xyz(unique_cations)


# In[11]:


df=pd.DataFrame.from_dict(unique_cation_dict, orient='index')
df=df.reset_index()
df=df.rename(columns={'index':'smiles'})
df.to_csv(x_file+'/'+x_file+'_rdkit_feature.csv',index=False)


# In[12]:


# print('dictionary creation finished!!!!')


# In[ ]:





# In[ ]:




