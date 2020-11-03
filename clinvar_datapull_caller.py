import json
import pickle
import os
from Bio import Entrez
# from dotenv import load_dotenv
from pathlib import Path
from clinvar_datapull import ClinVar_Datapull
from clinvar_parse import ClinVar_Parse
from clinvar_datapull_constants import DATAFILES_PATH
from drive_api import Drive_Api

# max number of retries after failed HTTP failures
Entrez.max_tries = 5
eutils_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
secrets_path = Path("~/secrets")

if __name__ == "__main__":
    # run pull and parse for cardio and cancer gene panels
    # and creates sheets from parse csvs
    data = ClinVar_Datapull(cardio=True, cancer=True)
    data.get_records()
    parse = ClinVar_Parse()
    parse.get_variation_data()
    # gene_pmids = parse.count_pmids()
    parse.convert_parse_to_csv()
    drive = Drive_Api()
    drive.create_sheets_from_csvs()
# TODO write argparse arguments
