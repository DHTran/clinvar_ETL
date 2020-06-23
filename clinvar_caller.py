import json
import pickle
import os
from Bio import Entrez
from collections import Counter
from datetime import date
# from dotenv import load_dotenv
from googleapiclient.discovery import build
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request
from pathlib import Path
from clinvar_datapull import ClinVar_Datapull
from clinvar_parse import ClinVar_Parse
from clinvar_datapull_constants import DATAFILES_PATH
from drive_api import Drive_Api
# max number of retries after failed HTTP failures
Entrez.max_tries = 5
eutils_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
secrets_path = Path("~/secrets")

def count_pmids():
    """takes list of dicts of Clinvar data where 
    [{'variation_id1': id1, 'pmid': [pmid1, pmid2],...},
     {'variation_id2': id2, 'pmid': []}]
    and counts occurrences of a given pmid number

    """
    gene_pmids = {}
    jsons_path = DATAFILES_PATH/'clinvar_parses'
    csvs_path = DATAFILES_PATH/'parse_csvs'
    for gene, parse in load_gene_json(jsons_path, load_parses=True):
        pmids = []
        for variation in parse:
            [pmids.append(x) for x in variation['pmids']]
        pmids_count = Counter(pmids).most_common()
        gene_pmids[gene] = pmids_count
    return gene_pmids

def load_gene_json(path, load_parses=False):
    """loads jsons for parsing
    """
    if load_parses:
        # load parses for conversion to other formats
        filter_files = '*_parse.json'
    else: 
        # load datapull jsons
        filter_files = '*.json'
    for file_ in Path(path).glob(filter_files):
        if load_parses:
            gene = file_.stem[:-6]
        else:
            gene = file_.stem
        if gene in parse.gene_list:  
            print(f"loading: {file_}")
            with open(file_) as f:
                data = json.load(f)
            yield gene, data
        else: 
            continue


data = ClinVar_Datapull(
    cardio=True, track_progress=True, overwrite=True)
data.get_records()
parse = ClinVar_Parse(overwrite=False)
gene_pmids = parse.count_pmids()
parse.convert_parse_to_csv()
drive = Drive_Api()
drive.create_sheets_from_csvs()