# clinvar_ETL

### Introduction 
Query NCBI ClinVar by gene(s) to pull all, associated non-CNV (by default)
records associated. Alternatvely pull records by variation id(s). 
- Save variation records as json data
- Read saved records
- parse ClinVar xmls (with Beautifulsoup library) to pull specific data (variant,
  submission, and classification data)
- Save parsed data as jsons or csvs
- Create Google sheets from csvs containing the ClinVar data

### Package structure

```.
├── Pipfile
├── Pipfile.lock
├── README.md
├── __init__.py
├── clinvar_ETL/
│   ├── __init__.py
│   ├── blocklist_pmids.csv
│   ├── clinvar_datapull.py
│   ├── clinvar_parse.py
│   ├── clinvar_utilities.py
│   ├── drive_api.py
│   ├── glaucoma_genes.csv
│   ├── main.py
│   ├── regex_patterns.py
│   ├── setup.py
│   └── tests/
│       └── test_datapull.py
├── clinvar_ETL_package.yml
├── poetry.lock
└── pyproject.toml
```

**Requires a datafiles folder to store jsons/csvs**  
~/datafiles/clinvar_ETL_datafiles/datapull_jsons/  
~/datafiles/clinvar_ETL_datafiles/parse_jsons/  
~/datafiles/clinvar_ETL_datafiles/parse_csvs/  


**Requires a .env file to  contains the following:**  
- Entrez.api_key = 'API_KEY' # to access 
- Entrez.email = 'your_email'
  - See https://www.ncbi.nlm.nih.gov/books/NBK25497/ for creation of an API key to facilitate access to NCBI databases
- A path to a secrets folder that loads Google Drive API credentials

### Usage

```
main.py provides a CLI for the package
-- main.py -o --gene 'MUTYH' --datapull:  queries ClinVar for MUTYH 
  variation records and save as MUTYH.json, -o: overwrite existing json file
-- main.py --datapull: query using default gene list (provided in 
  utilities/clinvar_utilities.py)
-- main.py --parse: loads datapull jsons and parses Variation records, 
  saves with 'gene_parse.json' filename pattern
-- main.py --validate: checks list of variation ids pulled using gene list
  versus gene.jsons stored in datapull folder
-- main.py --to_csv: converts parse_jsons to csv files (saved as 
  gene_parse_MM-DD-YYYY.csv
-- main.py --create_sheets: creates Google Sheets from csv files
-- main.py --update_sheets: update Google sheets formating to make it more
  readable
-- main.py --datapull --check:  Check arguments, don't run. 
```
