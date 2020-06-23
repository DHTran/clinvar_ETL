import csv
import json
import os
from Bio import Entrez
from pathlib import Path
from clinvar_ETL_constants import DATABASE, GENE_FILES
from clinvar_ETL_constants import DATAFILES_PATH, MULTIGENES
from dotenv import load_dotenv
load_dotenv()
Entrez.api_key = os.getenv('Entrez.api_key')
Entrez.email = os.getenv('Entrez.email')


class ClinVar_Datapull:
    """Pull ClinVar variation data for set of genes and returns selected
    data
    
    - uses a genes list from csv (cancer, cardio flags) or as a list 
      passed as an argument
    - pulls variation ids for a given gene using gene[Gene] queries
    with the optional single_gene flag
    - fetch_records pulls Variation page as text using variation ids
    - saves {'gene': [(id, record) tuples]} to json per gene

    -cancer, cardio, or both set to true loads gene names from the
    respective csvs with the defined gene lists

    -datafiles path is set to '~/DATAFILES_PATH' by default
    
    -test_flag = True limits records fetch to first 10, and returns 
    the id, records dict

    -overwrite=False (default) 
    prevents saving json of datapull if gene.json exists
    """
    
    def __init__(self, genes=None, test_flag=False, 
                 cancer=False, cardio=False, path=None,
                 overwrite=False, track_progress=False):
        if path:
            self.path = path
        else:
            self.path = (
                DATAFILES_PATH/'clinvar_datapull_datafiles/datapull_jsons')
        if genes:
            if isinstance(genes, list):
                self.genes = genes.upper()
            else:
                self.genes = genes.upper().split(",")
        elif cancer and cardio: 
            self.genes = self.create_gene_list(GENE_FILES, DATAFILES_PATH)
        elif cancer: 
            self.genes = self.create_gene_list(GENE_FILES[0], DATAFILES_PATH)
        elif cardio:
            self.genes = self.create_gene_list(GENE_FILES[1], DATAFILES_PATH)
        self.test_flag = test_flag 
        self.overwrite = overwrite
        self.track_progress = track_progress
        
    def __repr__(self):
        repr_string = (
            f"""genes {self.genes}, test_flag {self.test_flag}  
            path {self.path}, 
            overwrite is {self.overwrite}
            track_progress is {self.track_progress}"""
        )
        return repr_string
        
    def get_records(self):
        """gets records from ClinVar as
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
            gene_records = self.ids_by_gene(gene)
            if self.test_flag:
                return gene_records
            else:
                self.store_data_as_json(
                    gene_records, self.path, gene_json)
        
    def ids_by_gene(self, gene):
        """queries ClinVar by gene with gene[Gene]
        option: single_gene[prop] - filters to single gene results
        """
        if self.track_progress:
            print(f"in ids_by_gene for {gene}")
        gene_records = {}
        # some genes (e.g. RAD51D are listed as ≥2 genes)
        if gene not in MULTIGENES:
            terms = f'{gene}[Gene], single_gene[prop]'
        else:
            terms = f'{gene}[Gene]'
        print(terms)
        ids = self.get_ids(terms)
        if self.test_flag:
            ids = ids[0:10]
        id_records = self.fetch_clinvar_records(ids)
        gene_records[gene] = id_records
        # records is list of tuples: [(id_, record)...]  
        return gene_records

    def get_ids(self, terms):
        """gets database ids using search terms
        """
        if self.track_progress:
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
    def create_gene_list(cls, gene_files, path):
        """creates list of genes from csv files.  gene_files is
        list of paths to csv listing genes
        if test_genes=True, then uses test_genes list instead
        """
        genes = []
        if isinstance(gene_files, list):
            pass
        else: 
            gene_files = [gene_files]
        for file_ in gene_files:
            with open(path/file_) as f:
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
        """converts to json to store data locally
        """
        print(f"in store_data_as_json, storing {filename}")
        with open(path/filename, 'w') as f:
            json.dump(data, f)
    