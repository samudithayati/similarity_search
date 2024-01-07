# Create a Database

Create a database for molecules by generating properties using RDKIT. The rdkit is used to generate 17 properties in the python [script](https://github.com/samudithayati/similarity_search/blob/main/dataset_creation/dictionary_creation.py) for CSV files. 

## Usage

python3 [dictionary_creation.py](https://github.com/samudithayati/similarity_search/blob/main/dataset_creation/dictionary_creation.py) [001.csv](https://github.com/samudithayati/similarity_search/blob/main/dataset_creation/001.csv)
```
python3 dictionary_creation.py 001.csv
```

use python script [bash_create.py](https://github.com/samudithayati/similarity_search/blob/main/dataset_creation/bash_create.py) to create bash submission files for each CSV files. 

```
for x in {1..5};do sbatch $x.sh;done
```
