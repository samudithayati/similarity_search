# similarity_search

This project is based on finding similar molecules about a given molecule/molecules.
A ~18 million molecules database was compiled using Chembl and other datasets. The SMILES strings of these molecules were used in creating the database. Based on the SMILES string 200 embeddings were generated for each molecule using GROVER model.  These embeddings were used to find the similarity between two molecules. The cosine and L2 distances were used to rank the molecules in the database about the given molecule.   


## Usage
Use a csv file with all the smiles you need to search for example 
![result](https://github.com/samudithayati/similarity_search/example.csv)
