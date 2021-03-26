import csv
import json
import simplejson as sjson
import os
import pickle
import time
from bs4 import BeautifulSoup as bsoup
from Bio import Entrez
from pathlib import Path
from clinvar_ETL_constants import DATABASE, CANCER_GENE_SHEET
from clinvar_ETL_constants import ACMG59_GENE_SHEET, CARDIO_GENE_SHEET
from clinvar_ETL_constants import DATAFILES_PATH, MULTIGENES
from clinvar_ETL_constants import TESTS_PATH, TEST_RECORDS_IDS
from clinvar_ETL_constants import TEST_RECORDS_PATH
from dotenv import load_dotenv
load_dotenv()
Entrez.api_key = os.getenv('Entrez.api_key')
Entrez.email = os.getenv('Entrez.email')
# maxnumber of retries after failed HTTP requests
Entrez.max_tries = 5
secrets_path = Path("~/secrets")


class EncodeJsonMixin:
    def sjson_store(self, data, path):
        """store data locally with simplejson, assumes path is pathlib
        object
        """
        print(f"storing {path.name}")
        try:
            with open(path, 'w') as file_:
                sjson.dump(data, file_)
        except TypeError as e:
            print(f"TypeError saving data of type {type(data)}")
            print(e)

    def json_store(self, data, path):
        """store data locally with json. assumes path is pathlib
        object

        by default using json using compact separators option
        (','; ':')
        """
        print(f"storing {path.name}")
        try:
            with open(path, 'w') as file_:
                json.dump(data, file_, separators=(',', ':'))
        except TypeError as e:
            print(f"TypeError saving data of type {type(data)}")
            print(e)


class ClinVar_Datapull(EncodeJsonMixin):
    """Pull ClinVar variation data for set of genes and returns selected
    data

    Args:
    - gene_panel:  identifies gene panel to send to create_gene_list
        function to pull genes from csv(s) [e.g. CARDIO_GENE_LIST].
        If gene_panel = 'cardio' then pulls genes from 'cardio.csv'.
        Possible choices = 'cancer','cardio', and 'acmg59'.  Default is
        both cancer and cardio genes csvs if genes argument is none
    - genes:  list of genes if a custom list is desired.  default = None
    - test_flag (bool): True limits records fetch to first 10, and
        returns the id, records dict without saving, default is True
    - return_data (bool): if True returns data instead of saving
    - path (str) is set to '~/DATAFILES_PATH' by default
    - overwrite (bool): True allows overwriting of existing jsons,
        default is False

    main caller: get_records()
        pulls variation ids for a given gene using
        gene[Gene] queries with the optional [single_gene]
        property on the NCBI database
    """

    def __init__(self, gene_panel=None, genes=None,
                 return_data=False, path=None,
                 overwrite=False, test_flag=True):
        if path:
            self.path = path
        else:
            self.path = (
                DATAFILES_PATH/'datapull_jsons')
        self.gene_panel = gene_panel
        if genes is None:
            self.genes = self.create_gene_list(self.gene_panel)
        else:
            if isinstance(genes, list):
                self.genes = genes
            else:
                self.genes = [genes]
        self.return_data = return_data
        self.overwrite = overwrite
        self.test_flag = test_flag

    def __repr__(self):
        repr_string = (f"gene panel: {self.gene_panel} \
            \ngenes = {self.genes} \
            \nreturn_data: {self.return_data} \
            \npath: {self.path} \
            \ntest_flag: {self.test_flag} \
            \noverwrite: {self.overwrite}")
        return repr_string

    def get_records(self):
        """
        Iterates over self.gene (list) to query records from ClinVar as
        {'gene': [(id1, record1),(id2, record2)...]}

        Note a dict/json with all results is > 1 GB, saving each
        result per gene

        saves to self.path as gene.json.
        """
        for gene in self.genes:
            filename = f"{gene}.json"
            file_exists = self.file_exists(filename, self.path)
            if file_exists and not self.overwrite:
                print(f"skipping {filename}")
                continue
            gene_records = self.records_by_gene(gene)
            if self.return_data:
                return gene_records
            else:
                print(f"path = {self.path}")
                file_path = Path(self.path/filename)
                # store records for the given gene as json
                self.json_store(gene_records, file_path)

    def records_by_gene(self, gene):
        """
        Args:
            gene: name of gene to query ClinVar with 'gene[Gene]',
            'single_gene[prop]': filters to single gene results

        Return:
            gene_records: list of records for given gene
        """
        print(f"in records_by_gene for {gene}")
        # some genes (e.g. RAD51D are listed as ≥2 genes)
        if gene not in MULTIGENES:
            terms = f'{gene}[Gene], single_gene[prop]'
        else:
            terms = f'{gene}[Gene]'
        print(terms)
        ids = self.get_ids(terms)
        if self.test_flag:
            ids = ids[0:10]
        # ids records is list of tuples: [(id_, record)...]
        gene_records = self.fetch_clinvar_records(ids)
        return gene_records

    def get_ids(self, terms):
        """gets database ids using search terms
        """
        print("in get_ids")
        esearch = Entrez.read(Entrez.esearch(
            db=DATABASE,
            term=terms,
            retmax=100000,)
            )
        ids = esearch['IdList']
        count = esearch['Count']
        print(f"count {count}")
        return ids

    def fetch_clinvar_records(self, ids):
        """uses Entrez.efetch to pull Variation page
        data from ClinVar (rettype='vcv') as text

        https://biopython.org/docs/1.75/api/Bio.Entrez.html?highlight=efetch#Bio.Entrez.efetch
        https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch

        Args:
            ids = variation ids (numeric) for the variation pages
                  in ClinVar

        Return:
            records = list of records pulled with ids read with
            xml_reader method (default is txt)
        """
        records = []
        total = len(ids)
        for index, variation_id in enumerate(ids):
            print(f"fetching: id = {variation_id}, {index+1} of {total}")
            handle = self.fetch(
                db=DATABASE, id=variation_id, is_variationid=True,
                rettype='vcv', retmode='xml')
            record = self.xml_reader(handle)
            handle.close()
            if 'copy number' in record:
                print(f"copy number for {variation_id}")
                # avoids fetching CNV records
                continue
            records.append(record)
        return records

    def fetch(self, **kwargs):
        """method used for fetching ClinVar records

        ddefault = Entrez.efetch

        kwargs = (db=DATABASE, id=variation_id,
                is_variationid=True,
                rettype='vcv', retmode='xml')
        """
        time.sleep(0.1)
        return Entrez.efetch(**kwargs)

    def xml_reader(self, handle):
        """method to pull record from ClinVar with Biopython Entrez.fetch

        Args:
        - handle = data pulled from ClinVar (default = xml)
        Returns ClinVar record as txt
        """
        return handle.readlines()

    def filter_gene_list(self, user_input=False):
        """remove genes from genes list, filenames (gene.json) in the
        completed folder (option user_input)
        """
        genes_to_remove = []
        path = self.path
        if user_input:
            help_string = (
                "Enter genes to remove, separated by spaces or commas, \
                e.g. apc nbn brip1:  ")
            genes_to_remove = input(help_string)
            genes_to_remove = genes_to_remove.strip(",").upper().split(" ")
            print(f"removed {genes_to_remove}")
        else:
            for file in Path(path).glob('*.json'):
                gene = file.name.split(".")[0]
                genes_to_remove.append(gene)
        if isinstance(genes_to_remove, list):
            filtered_genes = [
                x for x in self.genes if x not in genes_to_remove]
            print(f"removed completed genes {genes_to_remove}")
            print(f"genes is now {filtered_genes}")
        self.genes = filtered_genes

    @classmethod
    def create_gene_list(cls, gene_panel):
        """creates list of genes from csv files. gene_files is
        list of paths to csv listing genes
        if test_genes=True, then uses test_genes list instead
        """
        genes = []
        gene_files = []
        if gene_panel == 'cancer':
            gene_files.append(CANCER_GENE_SHEET)
        elif gene_panel == 'cardio':
            gene_files.append(CARDIO_GENE_SHEET)
        elif gene_panel == 'acmg59':
            gene_files.append(ACMG59_GENE_SHEET)
        else:
            gene_files = [CANCER_GENE_SHEET, CARDIO_GENE_SHEET]
        # print(gene_files)
        for file_ in gene_files:
            with open(DATAFILES_PATH/file_) as f:
                reader = csv.reader(f)
                for row in reader:
                    gene = row[0].strip()
                    genes.append(gene)
        return genes

    def file_exists(self, filename, path):
        gene_json = Path(path/filename)
        print(f"checking for {gene_json}")
        if gene_json.is_file():
            print("file exists")
            return True
        else:
            return False

    @staticmethod
    def create_list_from_csvs(csvs, path):
        """creates list of items in first column of csvs. list_of_csvs is
        list of paths to the csv files
        """
        items = []
        for file in csvs:
            with open(path/file) as f:
                reader = csv.reader(f)
                for row in reader:
                    item = row[0].strip()
                    items.append(item)
        return items


class Test_records(ClinVar_Datapull):
    """
    ClinVar_Datapull subclass to create and store set of Clinvar
    records to test and validate.

    Args:
    test_ids (list of numbers) = list of variation ids to pull from
    ClinVar

    if test_ids = None uses TEST_RECORDS_IDS import from
    clinvar_ETL_constants.py

    """
    def __init__(self, test_ids=None, soup=False):
        ClinVar_Datapull.__init__(self, path=TESTS_PATH)
        if test_ids:
            self.test_ids = test_ids
        else:
            self.test_ids = TEST_RECORDS_IDS
        self.soup = soup

    def __repr__(self):
        repr_string = (f"""
            test_ids: {self.test_ids}
            soup flag: {self.soup}
        """)
        return repr_string

    def fetch_test_records(self, save=False):
        """
        method to fetch as set of test records using preset variation
        ids (TEST_RECORDS_IDS from clinvar_ETL_constants) or a
        specified set of ids
        """
        ids_records = None
        ids_records = self.fetch_clinvar_records(self.test_ids)
        if save:
            file_path = Path(self.path/'test.json')
            self.json_store(ids_records, file_path)
        else:
            return ids_records

    def soup_to_textfile(self, soup_records):
        """function to convert a list of soup objects
        to pickle files for local storage

        Args:
            soup_records = list of bsoup parsed xmls

        Return: list of paths to pickle files
        """
        file_paths = []
        for item in soup_records:
            variation_id = item.find("variationarchive").attrs['variationid']
            filename = f"{variation_id}_record"
            file_path = TEST_RECORDS_PATH/filename
            file_paths.append(file_path)
            with open(file_path, 'wb') as f:
                pickle.dump(item, f, protocol=pickle.HIGHEST_PROTOCOL)
        return file_paths

    def read_textfile_as_soup(self, file_paths):
        records = []
        for file_path in file_paths:
            print(file_path)
            with open(file_path, 'rb') as f:
                soup_record = pickle.load(f)
            records.append(soup_record)
        return records
