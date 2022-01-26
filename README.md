# bioassay-db
This repository contains the scripts to download Bioassays from PubChem and convert the results to a readable .csv file.
Check the associated [blogpost](https://medium.com/ersiliaio/downloading-pubchem-bioassays-made-easy-132ddd8e21c4) for more detailed explanations!

## Requirements
```
pandas
json
requests
```

## Installation
Create a conda environment (with Anaconda3 distribution, required packages are already installed) and clone this repository
```
conda create -n bioassays
conda activate bioassays
git clone https://github.com/ersilia-os/bioassay-db.git
```
For Windows users, please check [https://github.com/conda/conda/issues/8273] if you encounter errors with the SSL module

## Usage
Bioassays are downloaded as .json files (full information), and are later converted to .csv files with accompanying .txt description.
#### Download .json file
Add the desired AID_NUMBER (without any prefix), and the path to the destination folder. The download_json.py provides an example.
```
from src.pubchem import PubChemBioAssayRecord
c = PubChemBioAssayRecord(AID_NUMBER)
c.save_json("./examples")
```

#### Convert to csv
Add the folder path where the .json files are stored, as well as the assay to convert. This code snippet saves two different files, a .csv with all results and the assay description (including detailed explanations on each results column) in .txt format. The save_df.py provides an example.

Optionally, replace the function `get_all_results` by `get_outcome` to obtain a simplified version of the .csv results file containing solely the outcome of the bioassay.
```
from src.json2df import PubChemBioAssayJsonConverter
c = PubChemBioAssayJsonConverter("./examples", "PUBCHEM400.json")
df = c.get_all_results()
c.save_df(df, "./examples")
c.get_description("./examples")
```
