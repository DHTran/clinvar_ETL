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

```clinvar_ETL
├── main.py
├── Pipfile
├── Pipfile.lock
├── poetry.lock
├── pyproject.toml
├── README.md
├── datapull
│   └── clinvar_datapull.py
├── drive_api
│   └── drive_api.py
├── parse
│   └── clinvar_parse.py
├── tests
│   └── test_datapull.py
└── utilities
│     └── cclinvar_utilities.py
│     └── cglaucoma_genes.csv
│     └── regex_patterns.py
```

### Usage
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
