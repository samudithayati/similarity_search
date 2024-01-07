import pandas as pd
import glob
from tqdm import tqdm
#os.remove('tmp.csv')
       
for ii in tqdm(glob.glob('*.csv')):
    file = pd.read_csv(ii)
    #print("\nOriginal file:")
    #print(file)

    # adding header
    headerList = ['smiles']

    # converting data frame to csv
    file.to_csv(ii, header=headerList, index=False)
    
