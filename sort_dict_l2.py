import pickle
import pandas as pd
import sys

filename=sys.argv[1]
fname,ext=filename.split('.')
with open(filename, "rb") as ft1:   #Pickleloading
    rdkit_dict=pickle.load(ft1)

aa=list( rdkit_dict.keys())
name={}
for fn,ii in enumerate(aa):
    l2steve1=rdkit_dict[ii]['l2']
    new={k: v for k, v in sorted(l2steve1.items(), key=lambda item: item[1])}
    newd=pd.DataFrame(new.items(), columns=['molecule', 'l2'])
    alldf=newd[:100]
    alldf.loc[-1]=[ii,'reference']
    alldf.index = alldf.index + 1
    alldf = alldf.sort_index()
    alldf.to_csv(fname+'_l2_mol_'+str(fn)+'.csv',index=False)
    name[fn]=ii

with open(fname+'l2_naming.pickle','wb') as hp:
    pickle.dump(name, hp)