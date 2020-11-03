import csv
import simplejson as sjson
import os
from Bio import Entrez
from pathlib import Path
from clinvar_ETL_constants import DATABASE, CANCER_GENE_SHEET
from clinvar_ETL_constants import ACMG59_GENE_SHEET, CARDIO_GENE_SHEET
from clinvar_ETL_constants import DATAFILES_PATH, MULTIGENES
from clinvar_ETL_constants import TESTS_PATH, TEST_RECORDS_IDS
from dotenv import load_dotenv
load_dotenv()
Entrez.api_key = os.getenv('Entrez.api_key')
Entrez.email = os.getenv('Entrez.email')


class ClinVar_Datapull:
    """Pull ClinVar variation data for set of genes and returns selected
    data

    Args:
        gene_panel:  identifies gene panel to send to
            create_gene_list function to pull genes from csv(s)
            [e.g. CARDIO_GENE_LIST].  If gene_panel = 'cardio' then
            pulls genes from 'cardio.csv'.  Possible choice = 'cancer',
            'cardio', and 'acmg59'.  Default is both cancer and cardio
            genes csvs.

        genes:  list of genes if a custom list is desired.  default =
            None

        test_flag (bool): True limits records fetch to first 10, and returns
            the id, records dict without saving

        path (str) is set to '~/DATAFILES_PATH' by default

        overwrite (bool): True allows overwriting of existing jsons,
            default is False

    method get_records:
        pulls variation ids for a given gene using
            method ids_by_gene:
                gene[Gene] queries with the optional [single_gene]
                property on the NCBI database

    method fetch_records pulls Variation page as text using variation ids
    saves {'gene': [(id, record) tuples]} to json per gene
    """

    def __init__(self, gene_panel=None, genes=None, test_flag=False, path=None,
                 overwrite=False):
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
        self.test_flag = test_flag
        self.overwrite = overwrite

    def __repr__(self):
        repr_string = (f"""
            gene panel: {self.gene_panel}
            genes = {self.genes}
            test_flag: {self.test_flag}
            path: {self.path}
            overwrite: {self.overwrite}"""
        )
        return repr_string

    def get_records(self):
        """
        Args:
            None

        Iterates over self.gene (list) to query records from ClinVar as
        {'gene': [(id1, record1),(id2, record2)...]}

        Note a dict/json with all results is > 1 GB, saving each
        result per gene

        saves to self.path as gene.json.
        """
        for gene in self.genes:
            gene_json = f"{gene}.json"
            file_exists = self.file_exists(gene_json, self.path)
            if file_exists and not self.overwrite:
                print(f"skipping {gene_json}")
                continue
            gene_records = {}
            gene_records = self.records_by_gene(gene)
            if self.test_flag:
                return gene_records
            else:
                print(f"path = {self.path}")
                self.store_data_as_json(
                    gene_records, self.path, gene_json)

    def records_by_gene(self, gene):
        """
        Args:
            gene (list): list of genes to iterate over

        queries ClinVar by gene with gene[Gene]
        option: single_gene[prop] - filters to single gene results
        """
        print(f"in records_by_gene for {gene}")
        gene_records = {}
        # some genes (e.g. RAD51D are listed as ≥2 genes)
        if gene not in MULTIGENES:
            terms = f'{gene}[Gene], single_gene[prop]'
        else:
            terms = f'{gene}[Gene]'
        print(terms)
        ids = self.get_ids(terms)
        if self.test_flag:
            ids = ids[0:9]
        # ids records is list of tuples: [(id_, record)...]
        id_records = self.fetch_clinvar_records(ids)
        gene_records[gene] = id_records
        return gene_records

    def get_ids(self, terms):
        """gets database ids using search terms
        """
        print(f"in get_ids")
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
        """
        ids_records = []
        total = len(ids)
        for index, variation_id in enumerate(ids):
            id_record = {}
            print(f"fetching: id = {variation_id}, {index+1} of {total}")
            handle = Entrez.efetch(
                db=DATABASE, id=variation_id, is_variationid=True,
                rettype='vcv', retmode='xml')
            record = handle.readlines()
            if 'copy number' in record:
                print(f"copy number for {variation_id}")
                handle.close()
                continue
            id_record[variation_id] = record
            ids_records.append(id_record)
            handle.close()
        return ids_records

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
        print(gene_files)
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
            print(f"file exists")
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

    @staticmethod
    def store_data_as_json(data, path, filename):
        """converts to json to store data locally (with simplejson)
        """
        print(f"in store_data_as_json, storing {filename}")
        try:
            with open(path/filename, 'w') as f:
                sjson.dump(data, f)
        except TypeError as e:
            print(f"TypeError saving data of type {type(data)}")

class Test_records(ClinVar_Datapull):
    """
    ClinVar_datapull subclass to create and store set of Clinvar
    records to test and validate

    Default set:
        https://www.ncbi.nlm.nih.gov/clinvar/variation/55634/
        - NM_007294.4(BRCA1):c.5576C>G (p.Pro1859Arg)
        - several ref.
        - LB/Benign


        https://www.ncbi.nlm.nih.gov/clinvar/variation/55628/
        - NM_007294.4(BRCA1):c.5558dup (p.Tyr1853Ter)
        - several ref.
        - pathogenic

        https://www.ncbi.nlm.nih.gov/clinvar/variation/127898/
        - multi-gene (RAD51L3-RFFL also listed)
        - several ref.
        - VUS

        https://www.ncbi.nlm.nih.gov/clinvar/variation/464706/
        - NM_001128425.1(MUTYH):c.1626C>T (p.His5
        - silent


        https://www.ncbi.nlm.nih.gov/clinvar/variation/37003/
        - CDKN2A, 5-BP DUP, NT19
        - unmapped
    """
    def __init__(self, test_ids=None):
        ClinVar_Datapull.__init__(self, path=TESTS_PATH)
        if test_ids:
            self.test_ids = test_ids
        else:
            self.test_ids = TEST_RECORDS_IDS


    def __repr__(self):
        repr_string = (f"""
            genes: {self.genes}

            test_flag: {self.test_flag}
            gene_list: {self.gene_files}
            path: {self.path}
            overwrite: {self.overwrite}

            test_ids: {self.test_ids}
        """
        )
        return repr_string

    def create_test_records(self):
        """
        method to create set of test records calling
        ClinVar_Datapull.fetch_clinvar_records(self.test_ids)
        """
        ids_records = None
        ids_records = (
            ClinVar_Datapull.fetch_clinvar_records(self, self.test_ids)
            )
        return ids_records

    def create_save_test_records(self):
        """
        method to create and save test_records as jsons

        calls ClinVar_Datapull.fetch_clinvar_records(self.test_ids) and
        ClinVar_Datapull.store_data_as_json(
            test_records, path, 'test_records.json')
        """
        ids_records = None
        ids_records = self.create_test_records()
        ClinVar_Datapull.store_data_as_json(
            ids_records, self.path, 'test_records.json')