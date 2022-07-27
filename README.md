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

```clinvar_ETL/
├── README.md
├── clinvar_ETL
│   ├── __init__.py
│   ├── __main__.py
│   ├── clinvar_datapull.py
│   ├── clinvar_parse.py
│   ├── clinvar_utilities.py
│   ├── drive_api.py
│   ├── clinvar_ETL
│   ├── regex_patterns.py
│   ├── datafiles
│   │   ├── blocklist_pmids.csv
│   │   ├── datapull_jsons
│   │   ├── glaucoma_genes.csv
│   │   ├── missing_ids
│   │   ├── parse_csvs
│   │   ├── parse_jsons
│   │   └── test_records
│   └── tests
│       └── test_datapull.py
├── poetry.lock
```

**Requires a datafiles files and folders**
~/clinvar_ETL/clinvar_ETL/datafiles/datapull_jsons/  
~/clinvar_ETL/clinvar_ETL/datafiles/parse_jsons/  
~/clinvar_ETL/clinvar_ETL/datafiles/parse_csvs/  
~/clinvar_ETL/clinvar_ETL/datafiles/blocklist_pmids.csv  
~/clinvar_ETL/clinvar_ETL/datafiles/glaucoma_genes.csv  

**Requires a .env file to  contains the following:**
- Entrez.api_key = 'API_KEY'
- Entrez.email = 'your_email'
  - See https://www.ncbi.nlm.nih.gov/books/NBK25497/ for creation of an API key to facilitate access to NCBI databases
- A path to a secrets folder that loads Google Drive API credentials
  - 7/22/22: drive_api.py needs updating to work with recent Drive API updates

### Usage

```
/clinvar_ETL/clinvar_ETL/main.py provides a CLI for the package and can be
run directly or via python clinvar_ETL
-- clinvar_ETL -o --gene 'MUTYH' --datapull:  queries ClinVar for MUTYH
  variation records and save as MUTYH.json, -o: overwrite existing json file
-- clinvar_ETL --datapull: query using default gene list (provided in
  utilities/clinvar_utilities.py)
-- clinvar_ETL --parse: loads datapull jsons and parses Variation records,
  saves with 'gene_parse.json' filename pattern
-- clinvar_ETL --validate: checks list of variation ids pulled using gene list
  versus gene.jsons stored in datapull folder
-- clinvar_ETL --to_csv: converts parse_jsons to csv files (saved as
  gene_parse_MM-DD-YYYY.csv
-- clinvar_ETL --create_sheets: creates Google Sheets from csv files
-- clinvar_ETL --update_sheets: update Google sheets formating to make it more
  readable
-- clinvar_ETL --datapull --check:  Check arguments, don't run.
```
