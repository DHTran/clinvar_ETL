"""Pull data from ClinVar (variation records) and stores as jsons.
"""

import csv
import json
# required to avoid "Couldn't find tree builder with the features..." error
# import lxml
import simplejson as sjson
import os
import pickle
import time
from bs4 import BeautifulSoup as bsoup
from Bio import Entrez
from pathlib import Path
from clinvar_ETL.clinvar_utilities import DATABASE, GLAUCOMA_GENES
from clinvar_ETL.clinvar_utilities import PACKAGE_DATAFILES, DATAPULL_JSONS_PATH
from clinvar_ETL.clinvar_utilities import TEST_RECORDS_PATH
from dotenv import load_dotenv
# from typing import Union
from urllib.error import HTTPError

load_dotenv()
Entrez.api_key = os.getenv('Entrez.api_key')
Entrez.email = os.getenv('Entrez.email')
# maxnumber of retries after failed HTTP requests
Entrez.max_tries = 5


class Decode_Jsons_Mixin:
    """Mixin to decode json files (ClinVar datapulls and parses).

    Note this Mixin is a duplicate of Mixin in clinvar_parse.py
    """

    def decode(self, file):
        """Retrieve json-stored data with simplejson."""
        try:
            data = sjson.load(file)
        except TypeError as e:
            print(e)
            print("Trying json instead")
            data = json.load(file)
        except ValueError as e:
            print(e)
            print("Trying json instead")
            data = json.load(file)
        return data

    def decodes(self, file):
        """Retrieve json-stored data with simplejson."""
        try:
            data = json.loads(file)
        except TypeError:
            print("Trying sjson instead")
            data = sjson.loads(file)
        return data

    def load_gene_json(self, path, genes=None, file_filter=None):
        """Load files based on file_filter and yield gene: data dicts
        (assumes json encodes a dict with {gene: clinvar_data}.

        Arg:
            file_filter -- str to filter directory with glob
            genes -- list of genes to check before attempting to decode
            path -- default is DATAPULLS_PATH
        """
        for filepath in Path(path).glob(file_filter):
            print(f"loading {filepath}")
            if 'parse' in filepath.stem:
                # extracts gene name from "gene_parse.json"
                gene = filepath.stem[:-6]
            else:
                # extracts gene name from "gene.json"
                gene = filepath.stem
            if gene in genes:
                # checks for gene.json with gene in panel list
                with open(filepath) as file:
                    if file:
                        try:
                            clinvar_data = self.decode(file)
                        except AttributeError:
                            clinvar_data = self.decodes(file)
                yield gene, clinvar_data
            else:
                continue


class EncodeJsonMixin:
    """Mixin to provide json/sjson dump methods."""

    def store(self, data, file):
        """Store data locally with sjson then json if TypeError occurs
        (ClinVar xml can fail).

        Args:
            data -- data to encode
            path -- path to store file (checks if Pathlib Path)
        """
        try:
            sjson.dump(data, file)
        except TypeError:
            json.dump(data, file)


class Fetch_Mixin:
    """Mixin to query ClinVar by variation ids and return xml records using the
    BioPython Entrex.fetch method to query the an NCBI database.
    """

    def fetch_clinvar_records(self, ids: list[int]):
        """Use Entrez.efetch to pull Variation page data from ClinVar
        (rettype='vcv') as text.

        https://biopython.org/docs/1.75/api/Bio.Entrez.html?highlight=efetch#Bio.Entrez.efetch
        https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch

        Args:
            ids = list of variation ids [should be list of integer(s)
            but may be string(s)] for the variation pages in ClinVar

        Return:
            records = list of records pulled with ids read with
            xml_reader method (default is txt)
        """
        records = []
        # convert comma separate string to list
        if isinstance(ids, str):
            ids = list(ids.split(','))
        # if ids is a list of integers
        if isinstance(ids, list):
            if all(isinstance(id, int) for id in ids):
                pass
            else:
                if not all(id.isnumeric() for id in ids):
                    raise Exception("ids contains non-digits")
                ids = [int(id) for id in ids]
        # note ClinVar truncates decimals apparently
        elif isinstance(ids, set) or isinstance(ids, tuple):
            ids = list(ids)
        else:
            # if ids is a single integer,
            ids = [ids]
        total = len(ids)
        for index, variation_id in enumerate(ids):
            print(f"fetching: id = {variation_id}, {index+1} of {total}")
            record = self.fetch(
                db=DATABASE, id=variation_id, is_variationid=True,
                rettype='vcv', retmode='xml')
            if record is None:
                record = ("warning empty record", variation_id)
            records.append(record)
        return records

    def fetch(self, **kwargs):
        """Pull ClinVar records using Biopython's efetch.

        kwargs -- (db=DATABASE, id=variation_id, is_variationid=True,
                   rettype='vcv', retmode='xml')
        """
        retries = 1
        while retries < 6:
            try:
                handle = Entrez.efetch(**kwargs)
                record = handle.readlines()
                handle.close()
                return record
            except HTTPError as e:
                print(f"HTTPerror {e}")
                print(f"{kwargs}")
                time.sleep(0.2)
                retries += 1
        variation_id = kwargs['id']
        return ('fetch error', variation_id)


class ClinVar_Datapull(EncodeJsonMixin, Fetch_Mixin):
    """Pull ClinVar variation data for set of genes and returns selected
    data.

    Args:
      gene_panel -- identifies gene panel to send to create_gene_list
        function to pull genes from csv(s) [e.g. CARDIO_GENE_LIST].
        If gene_panel = 'cardio' then pulls genes from 'cardio.csv'.
        Possible choices = 'cancer','cardio', 'acmg59', 'new_panel'.
        Default is both cancer and cardio genes csvs if genes argument is none
      genes --  list of genes if a custom list is desired.  default = None
      test_flag (bool) -- if True limits records fetch to first 10, and
        returns the id, records dict without saving, default is False
      return_data (bool) -- if True returns data instead of saving
      path (str) -- path to datafiles, is set to 'PACKAGE_DATAFILES' by default
      overwrite (bool) -- True allows overwriting of existing jsons,
        default is False

    main caller is get_records()
        pulls variation ids for a given gene using
        gene[Gene] queries with the optional [single_gene]
        property on the NCBI database
    """

    def __init__(self, gene_panel=None, genes=None,
                 return_data=False, path=None,
                 overwrite=False, test_flag=False):
        if path:
            self.path = path
        else:
            self.path = DATAPULL_JSONS_PATH
        self.gene_panel = gene_panel
        # uses gene_panel to build list of genes
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
        self.ids_pulled = []

    def __repr__(self):
        repr_string = (
            f"gene panel: {self.gene_panel} \n"
            f"genes = {self.genes} \n"
            f"return_data: {self.return_data} \n"
            f"path: {self.path} \n"
            f"test_flag: {self.test_flag} \n"
            f"overwrite: {self.overwrite} \n")
        return repr_string

    def get_records(self):
        """Iterate over self.gene (list) to query records from ClinVar as
        {'gene': [(id1, record1),(id2, record2)...]}

        Note a dict/json with all results is > 1 GB, saving each
        result per gene

        Saves to self.path as gene.json.
        """
        for gene in self.genes:
            if self.test_flag:
                filename = f"test_{gene}.json"
            else:
                filename = f"{gene}.json"
            file_exists = self.check_file_exists(filename, self.path)
            if file_exists and not self.overwrite and not self.return_data:
                print(f"skipping {filename}")
                continue
            gene_records = self.get_records_per_gene(gene)
            if gene_records:
                if self.return_data:
                    return gene_records
                else:
                    file_path = Path(self.path/filename)
                    # store records for the given gene as json
                    with open(file_path, 'w') as file:
                        print(f"storing {file_path}")
                        self.store(gene_records, file)

    def get_records_per_gene(self, gene):
        """Assemble search terms for passing to fetch_clinvar_records
        method.

        Returns list of records for a given gene.

        Args:
            gene -- name of gene to query ClinVar with 'gene[Gene]',
            'single_gene[prop]': filters to single gene results

        Return:
            gene_records -- list of records for given gene
        """
        print(f"records_by_gene for {gene}")
        terms = self.get_query_terms(gene)
        print(terms)
        gene_ids = self.get_ids_per_gene(terms)
        if gene_ids:
            if self.test_flag:
                gene_ids = gene_ids[0:10]
            gene_records = self.fetch_clinvar_records(gene_ids)
            return gene_records
        elif gene_ids is None:
            return None

    @staticmethod
    def get_query_terms(gene):
        "Pull ClinVar ids for a given gene."

        # filters out CNVs
        not_filters = [' NOT "copy number gain"', ' NOT "copy number loss"']
        terms = f'{gene}[Gene]'
        for not_tag in not_filters:
            terms += f'{not_tag}[Type of variation]'
        return terms

    @staticmethod
    def get_ids_per_gene(terms, retmax=100000, retstart=1):
        """Retrieve NCBI database (db=DATABASE) ids using search terms.

        Note retmax = 100000 (max records returned)
        current ClinVar limit is 100,000 per query but ClinVar records
        are >100,000. retstart default 1.
        """
        ids = []
        total_count = int(0)
        # may need to loop several times since query is capped at retmax
        while len(ids) <= total_count:
            # print(f"retstart: {retstart}")
            esearch = Entrez.read(
                Entrez.esearch(db=DATABASE, term=terms,
                               retmax=retmax, retstart=retstart))
            # esearch['Count'] is always total number
            total_count = int(esearch.get('Count'))
            if total_count == 0:
                # Query results return no results e.g. all CNVs for a gene
                print(f"{terms} returned {total_count} ids")
                return None
            query_ids = esearch['IdList']
            ids.extend(query_ids)
            retstart += retmax
            if len(ids)+1 == total_count:
                return ids
        return ids

    def filter_gene_list(self, user_input=False):
        """Remove genes from genes list, filenames (gene.json) in the
        completed folder (option user_input).
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
        """Create list of genes from csv files.

        Arg:
            gene_panel: if no gene_panel is set then gene_file = GLAUCOMA_GENES
            GLAUCOMA_GENES points to 'glaucoma_genes.csv' containing genes list

            # As of 9/2022, glaucoma is the only panel
        """
        genes = []
        if gene_panel == 'glaucoma':
            gene_file = GLAUCOMA_GENES
        elif gene_panel == 'TBD':
            # for future panels
            gene_file = None
        else:
            gene_file = GLAUCOMA_GENES
        with open(PACKAGE_DATAFILES/gene_file) as f:
            reader = csv.reader(f)
            for row in reader:
                gene = row[0].strip()
                genes.append(gene)
        return genes

    def fetch_test_records(self, save=False):
        """Fetch as set of test records using preset variationids
        (TEST_RECORDS_IDS from clinvar_utilities) or a specified set of ids.
        """
        test_ids = TEST_RECORDS_IDS
        ids_records = self.fetch_clinvar_records(test_ids)
        if save:
            file_path = Path(self.path/'test.json')
            self.store(ids_records, file_path)
        else:
            return ids_records

    @staticmethod
    def check_file_exists(filename, path):
        """Check if file exists with is_file()"""
        file_path = Path(path/filename)
        print(f"checking for {file_path}")
        if file_path.is_file():
            print("file exists")
            return True
        else:
            return False

    @staticmethod
    def create_list_from_csvs(csvs, path):
        """Create list of items in first column of csvs. list_of_csvs is
        list of paths to the csv files.
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
    def soup_to_textfile(self, soup_records):
        """Convert a list of soup objects to pickle files for local storage.

        Args:
            soup_records -- list of bsoup parsed xmls

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

    @staticmethod
    def read_pickled_file(file_paths):
        """Open pickled data."""
        records = []
        for file_path in file_paths:
            print(file_path)
            with open(file_path, 'rb') as f:
                soup_record = pickle.load(f)
            records.append(soup_record)
        return records


class Validate_Datapull(Fetch_Mixin, EncodeJsonMixin, Decode_Jsons_Mixin):
    """Compares Variation ids in saved jsons versus ids from [Gene] queries
    of the ClinVar database.

    Args:
        gene_panel --  identifies gene panel to send to create_gene_list
            function to pull genes from csv(s)
        genes --  list of genes if a custom list is desired.  default = None
            path: path to folder holding jsons default = DATAPULL_JSONS_PATH
        save_json (bool) -- bool, if True save missing records to json,
            default = False
        path (Path object) -- path to folder holding datapull jsons (default)
    """

    def __init__(self, genes=None, gene_panel=None, save_json=False, path=None):
        self.gene_panel = gene_panel
        if genes is None:
            self.genes = ClinVar_Datapull.create_gene_list(self.gene_panel)
        else:
            if isinstance(genes, list):
                self.genes = genes
            else:
                self.genes = [genes]
        if path is None:
            self.path = DATAPULL_JSONS_PATH
        self.save_json = save_json

    def __repr__(self):
        repr_string = (f"""
        genes = {self.genes},
        path = {self.path},
        gene_panel = {self.gene_panel}
        save_json = {self.save_json}
        """)
        return repr_string

    def check_and_update(self):
        """Main caller for Validate_Datapull class.

        Get ClinVar Variation ids missing from datapull and updates the gene
        jsons with missing records

        Note datapull_only_ids (2nd item returned) are not used in this method
        """
        missing_query_ids, _ = self.compare_ids()
        if missing_query_ids is None:
            print("No missing ids")
            return
        missing_query_ids = [int(id) for id in missing_query_ids]
        print(f"missing_query_ids: {len(missing_query_ids)}")
        user_check = input("Do you wish to continue the update?")
        if user_check.lower() == 'y':
            pass
        else:
            return
        missing_records = self.fetch_clinvar_records(missing_query_ids)
        if self.save_json:
            missing_records_name = Path(self.path/'missing_records.json')
            missing_ids_name = Path(self.path/'missing_ids.json')
            with open(missing_records_name, 'w') as file:
                sjson.dump(missing_records, file)
            with open(missing_ids_name, 'w') as file:
                sjson.dump(missing_query_ids, file)
        self.add_missing_records(missing_records)

    def load_datapull_json(self):
        """Retrieve ClinVar records stored as jsons.

        ClinVar records are in the format: {'gene': [list of records]}
        (note: parsed from xml using Biopython)
        """
        for filepath in Path(self.path).glob('*.json'):
            # print(f"loading {filepath}")
            gene = filepath.stem
            if gene in self.genes:
                with open(filepath) as file:
                    if file:
                        clinvar_data = sjson.load(file)
                yield gene, clinvar_data
            else:
                continue

    def ids_from_datapull_jsons(self):
        """Read datapull jsons to retrive the ClinVar ids that were successfully
        pulled.

        Args:
            path --  path to folder holding the datapull jsons
        """
        datapull_ids = {}
        for gene, records in self.load_datapull_json():
            ids = []
            for record in records:
                record = bsoup(str(record), features='lxml')
                variation_id = record.find(
                    "variationarchive").attrs['variationid']
                ids.append(variation_id)
            datapull_ids[gene] = ids
        return datapull_ids

    def pull_clinvar_query_ids(self):
        """Query ClinVar for all variation ids associated with genes in
        self.genes.

        Returns dict where: {'gene1': [id1, id2...], 'gene2': [id1, id2, ...]}
        """
        all_ids = {}
        for gene in self.genes:
            terms = ClinVar_Datapull.get_query_terms(gene)
            variation_ids = ClinVar_Datapull.get_ids_per_gene(terms)
            if variation_ids:
                all_ids[gene] = variation_ids
            elif variation_ids is None:
                pass
        return all_ids

    def compare_ids(self):
        """Validate datapull by matching pulled ids versus initial id query.

        Looks for ClinVar variation ids that were no successully datapulled
        (missing_query_ids) and ids that are present in datapull jsons but not
        in the list of ids associated with the genes query (datapull_only_ids).

        missing_query_ids are likely from Http errors or recently added
        ClinVar pages

        datapull_only_ids are likely ClinVar pages that were recently deleted

        Returns Sets of variation ids
        """
        def values_as_set(items):
            values_as_set = set()
            print(items.keys())
            for value in items.values():
                values_as_set.update(value)
            return values_as_set

        # ids returned as dicts of {'gene': [id1, id2...]}
        query_ids = self.pull_clinvar_query_ids()
        assert query_ids is not None, "query ids returned None"
        datapull_ids = self.ids_from_datapull_jsons()
        assert datapull_ids is not None, "datapull ids returned None"
        print(f"datapull query id count {len(datapull_ids)}")
        print(f"query ids count {len(query_ids)}")
        query_ids_set = values_as_set(query_ids)
        datapull_ids_set = values_as_set(datapull_ids)
        if query_ids_set == datapull_ids_set:
            print("datapull ids match query ids")
        elif len(query_ids_set) > len(datapull_ids_set):
            print("more query ids than datapull ids")
        elif len(query_ids_set) < len(datapull_ids_set):
            print("less query ids than datapull ids")
        missing_query_ids = query_ids_set.difference(datapull_ids_set)
        datapull_only_ids = datapull_ids_set.difference(query_ids_set)
        return missing_query_ids, datapull_only_ids

    def add_missing_records(self, records):
        """Identify missing records (difference of original gene query ids
        versus ids obtained after datapull).

        Adds missing records to appropriate gene.json
        """

        def get_genes_from_records():
            """Create dict where keys are genes from the missing data and
            the values are the missing ClinVar record data (as list of records)
            """
            missing_data_by_genes = {}
            # account for single record instead of a nested list of records
            nonlocal records
            if all(isinstance(x, list) for x in records):
                pass
            else:
                records = [records]
            for record in records:
                soup_record = bsoup(str(record), features='lxml')
                gene = soup_record.find("gene").attrs['symbol']
                if gene not in missing_data_by_genes:
                    missing_data_by_genes[gene] = [record]
                else:
                    missing_data_by_genes[gene].extend([record])
            return missing_data_by_genes

        user_check = input("This will overwrite and update the gene.json files. \
                            Do you wish to continue? y/n")
        if user_check.lower() == 'y':
            pass
        else:
            return
        missing_data_by_genes = get_genes_from_records()
        if self.save_json:
            save_path = Path(DATAPULL_JSONS_PATH/'missing_data_by_genes.json')
            with open(save_path, 'w') as file:
                sjson.dump(missing_data_by_genes, file)
        for gene, missing_data in missing_data_by_genes.items():
            file_name = f"{gene}.json"
            filepath = Path(DATAPULL_JSONS_PATH/file_name)
            if filepath.is_file():
                with open(filepath) as file:
                    print(f"opening {filepath}")
                    clinvar_data = sjson.load(file)
                clinvar_data.extend(missing_data)
                with open(filepath, 'w') as file:
                    print(f"updating {filepath}")
                    sjson.dump(clinvar_data, file)