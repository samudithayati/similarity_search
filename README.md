# similarity_search

This project is based on finding similar molecules about a given molecule/molecules.
A ~18 million molecules database was compiled using Chembl and other datasets. The SMILES strings of these molecules were used in creating the database. Based on the SMILES string 200 embeddings were generated for each molecule using GROVER model.  These embeddings were used to find the similarity between two molecules. The cosine and L2 distances were used to rank the molecules in the database about the given molecule.   


## Usage
Use a [CSV](https://github.com/samudithayati/similarity_search/blob/main/example.csv) file with all the smiles you need to use in similarity search.
If you are using a SLURM in an HPC use the submission script to submit the jobs by changing the appropriate keywords and pointers.
Use the Python script to perform the distance calculation for each molecule in  the database.

The calculated distance will structured in a dictionary format. The dictionary key will be the SMILES strings. Two sorting scripts were used to sort the [l2](https://github.com/samudithayati/similarity_search/blob/main/sort_dict_l2.py) and [cosine](https://github.com/samudithayati/similarity_search/blob/main/sort_dict_cosine.py) distances. These scripts will output the similarity search result for given distance matrices in a CSV file for each molecule.  
