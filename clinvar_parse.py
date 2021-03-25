import csv
import json
import simplejson as sjson
import re
import pandas as pd
from bs4 import BeautifulSoup as bsoup
from datetime import date
from collections import Counter, defaultdict
from pathlib import Path
from clinvar_ETL_constants import DATAPULL_JSONS_PATH, BLOCKLIST_PMIDS
from clinvar_ETL_constants import PARSE_JSONS_PATH
from clinvar_datapull import ClinVar_Datapull as clinvar_dp


class ClinVar_Parse_Caller:
    """main caller for parsing ClinVar records from datapull jsons

    """

    def parse_caller(self):
        pass


class ClinVar_Parse_Record:
    """
    class to parse a ClinVar Variation record and return selected
    information

    Arg:
        record:  a ClinVar Variation record (as Beautiful Soup4 object)
    """

    def __init__(self, record):
        self.record = record

    def __repr__(self):
        record_header = self.record.find('variationarchive').attrs
        repr_string = f"record: {str(record_header)}"
        return repr_string

    def get_record_data(self):
        """
        caller to parse out the specific data for a given ClinVar
        Variation Record

        variant-level data
        - variation_id, datelastupdated, chrom, start, stop, variant
        name, ref, alt,

        submission data (per submission)
        - submitter_name, classification, inheritance, condition,
        pmid(s), comment
        """
        variant_data = self.get_variant_data()
        variant_data['submission_data'] = self.get_assertion_data()
        variation_summary_data = self.variation_summary_data(
            variant_data)
        return variation_summary_data

    def get_variant_data(self):
        """pull name, variation_id, start, stop, ref, alt, chrom, assembly
        from clinvar record, returns dict.
        """
        results_dict = {}
        allele_data = ['start', 'stop', 'referenceallelevcf',
                       'alternateallelevcf', 'chr', 'assembly']
        variation_id = (
            self.record.find(
                "variationarchive").attrs['variationid'])
        date_updated = (
            self.record.find(
                "variationarchive").attrs['datelastupdated'])
        num_submitters = (
            self.record.find(
                "variationarchive").attrs['numberofsubmissions'])
        name = self.record.find("simpleallele").find("name").string
        sequencelocation_tags = self.record.find_all(
            "sequencelocation", assembly="GRCh37")
        if len(sequencelocation_tags) > 1:
            variant_data = self.record.find_all(
                "sequencelocation", assembly="GRCh37")[1].attrs
        else:
            # when VCF-stype data is not present
            variant_data = self.record.find_all(
                "sequencelocation", assembly="GRCh37")[0].attrs
        for item in allele_data:
            results_dict[item] = variant_data.get(item)
        results_dict['variation_id'] = variation_id
        results_dict['name'] = name
        results_dict['ref'] = results_dict.pop('referenceallelevcf')
        results_dict['alt'] = results_dict.pop('alternateallelevcf')
        results_dict['date_updated'] = date_updated
        results_dict['num_submitters'] = num_submitters
        # rename chr to chrom to avoid chr() function name
        results_dict['chrom'] = results_dict.pop('chr')
        return results_dict

    def get_assertion_data(self):
        """parses clinicalassertion node to get submission data

        returns list of submission data
        """
        clinical_assertions = self.record.select("clinicalassertion")

        def get_submitter_data(assertions):
            for assertion in assertions:
                submitter_data = {}
                pmids = []
                scv_number = (
                    assertion.find("clinvaraccession").
                    attrs['accession'])
                submitter_data['datelastupdated'] = (
                    assertion.attrs['datelastupdated'])
                submitter_data['submittername'] = (
                    assertion.find("clinvaraccession").
                    attrs['submittername'])
                submitter_data['classification'] = (
                    assertion.find("description").string.strip())
                submitter_data['inheritance'] = (
                    assertion.find("origin").string.strip())
                if assertion.find_all("citation"):
                    for item in assertion.find_all("citation"):
                        if item.id:
                            pmid_num = item.id.string.strip()
                            pmids.append(pmid_num)
                submitter_data['pmids'] = pmids
                if assertion.find('comment'):
                    submitter_data['comment'] = (
                        assertion.find('comment').string.strip())
                if assertion.find("traitset").find('elementvalue'):
                    disease = (
                        assertion.find("traitset").
                        find('elementvalue').string)
                else:
                    disease = assertion.find("xref").attrs
                if isinstance(disease, dict):
                    submitter_data['condition'] = disease
                elif isinstance(disease, str):
                    submitter_data['condition'] = disease.strip()
                else:
                    raise Exception(
                            f"disease is not a dict or str: {disease}")
                yield submitter_data
        assertion_results = (
            list(get_submitter_data(clinical_assertions)))
        assert len(assertion_results) == len(clinical_assertions), \
            f"""
            mismatch between # of submitter_data dicts
            {len(assertion_results)} and # of assertions
            {len(clinical_assertions)}"""
        return assertion_results

    def variation_summary_data(self, data):
        """method to summarize the parsed clinvar data in a user
        friendly format

        Args: data = variant_data dict (from get_assertion_data()
        method)
        """
        def classification_count(classification_data):
            """method to get number of classifications per
            classification"""
            classification_count = Counter(
                submitter['classification'] for submitter
                in classification_data)
            # remove python chars
            classification_count = ClinVar_Parse_Record.clean_string(
                dict(classification_count))
            return classification_count

        def group_classifications(data):
            """method to group submitters by classification
                and save as dict:
                {'classification1': [list of submitters],
                'classification2: [list of submitters]}
            """
            def group_submitters(submission_data):
                """reorganize submission data from
                {'submittername': 'classification'}

                to {'classification': ['submitter1', 'submitter2'...]}

                Arg: submission_data = variant_data dict
                """
                classification_submitters = defaultdict(list)
                for submission in submission_data:
                    submitter = submission.get('submittername')
                    classification = submission.get('classification')
                    classification_submitters[classification].append(submitter)
                return classification_submitters

            def clean_classifications(classifications):
                """
                Remove python characters and group by classification

                Args:
                    classifications (list of dicts):
                    [{'submittername': classification}, ...]

                Return:
                    string = 'classifications': submitter1, submitter2...
                """
                result_string = ''
                for classification, submitters in classifications.items():
                    class_submitters = ''
                    classification = classification + ": "
                    submitters = (
                        ClinVar_Parse_Record.clean_string(submitters))
                    class_submitters = (
                        "*  " + classification + submitters + ",  ")
                    result_string = result_string + class_submitters
                return result_string

            classification_submitters = dict(
                group_submitters(data['submission_data']))
            cleaned_classifications = (
                clean_classifications(classification_submitters))
            return cleaned_classifications

        def get_submitter_comments(submission_data):
            """get comment data from submissions
            """
            comment_dict = {}
            for submission in submission_data:
                submitter = submission.get('submittername')
                comment = submission.get('comment', 'None')
                if submitter not in comment_dict.keys():
                    comment_dict[submitter] = comment
                else:
                    comment_dict[submitter] += f", {comment}"
            return comment_dict

        def get_pmids(submission_data):
            """get pmids from submissions and flag (as bold) pmids
            on blocklist
            """
            pmids = []
            for submission in submission_data:
                for pmid in submission['pmids']:
                    if pmid not in pmids:
                        pmids.append(pmid)
                    else:
                        pass
            pmid_string = ''
            for pmid in pmids:
                if pmid in BLOCKLIST_PMIDS:
                    pmid = f"*{pmid}*"
                pmid = pmid + ", "
                pmid_string += pmid
            pmid_string = pmid_string[:-2]
            return pmid_string

        summary_dict = {}
        summary_dict['classification_count'] = (
            classification_count(data['submission_data']))
        summary_dict['classifications'] = group_classifications(data)
        summary_dict['comments'] = get_submitter_comments(
            data['submission_data'])
        summary_dict['pmids'] = get_pmids(data['submission_data'])
        keys_to_add = [
            'variation_id', 'name', 'date_updated',
            'num_submitters']
        for key in keys_to_add:
            summary_dict[key] = data.get(key)
        return summary_dict

    @staticmethod
    def clean_string(string):
        """simple method to remove data-specific characters
        for better user presentation
        """
        PYTHON_CHARS = re.compile("[\{\}\'\[\]]")

        if isinstance(string, str):
            pass
        else:
            string = str(string)
        return re.sub(PYTHON_CHARS, '', string)


class DecodeJsonMixin:
    """Mixin class to decode a json with option to
    use json or sjson library and choice of deserializer
    method

    Args:
        file_ = file object
        decode_type: 'load' or 'loads' default 'load' (eg. json.load)
    """
    def sjson_decode(self, file_, decode_type='load'):
        """retrieve json-stored data with simplejson
        """
        print(f"loading {file_}")
        if decode_type == 'load':
            data = sjson.load(file_)
        return data

    def json_decode(self, file_, decode_type='load'):
        """retrieve json-stored data. default json.load.
        Enter decode_type='loads' for json.loads
        """
        print(f"loading {file_}")
        if decode_type == 'load':
            data = json.load(file_)
        elif decode_type == 'loads':
            data = json.loads(file_)
        return data


class Read_jsons(DecodeJsonMixin):
    """read variation data pulled from ClinVar, stored as json

    Args:
        path: path to save parse jsons. Default is DATAPULL_JSONS_PATH
            from clinvar_ETL_constants.py

        return_data(bool): if True returns record(s) and count instead
            of saving to json. Default = False
        overwrite(bool): if True saves over existing file,
            Default = False
        genes(list): list of genes, default: pulls
            list of cancer and cardio genes (using imported
            ClinVar_Datapull class)
        suppress_out(bool): if True prints suppresses most print
            statements that log progress
        decoder: json package to use. default 'sjson' (alternative is
            json)
        decode_type: decode method from json package. default is 'load'
    """
    def __init__(self, path=None, overwrite=False, return_data=False,
                 genes=None, suppress_output=True, decoder='sjson',
                 decode_type='loadk'):
        if path:
            self.path = path
        else:
            self.path = DATAPULL_JSONS_PATH
        self.overwrite = overwrite
        self.return_data = return_data
        if genes is None:
            # ClinVar_Datapull by default creates genes list from
            # cancer and cardio lists
            genes = clinvar_dp().genes
            self.gene_list = genes
        else:
            if not isinstance(genes, list):
                genes = [genes]
            self.gene_list = genes
        self.suppress_output = suppress_output
        self.decoder = decoder
        self.decode_type = decode_type

    def __repr__(self):
        repr_string = (
            f"path {self.path}\n\
            return_data is {self.return_data}\n\
            overwrite is {self.overwrite}\n\
            gene_list = {self.gene_list}\n\
            blocklist_pmids = {BLOCKLIST_PMIDS}\n\
            decoder = {self.decoder}\n\
            decode_type = {self.decode_type}\n"
        )
        return repr_string

    def get_datapulls(self):
        """read json files storing variation records by gene from
        DATAPULL_JSONS_PATH (default), hard-coded filter '*.json'

        Return:  TODO
        """
        gene_records = {}
        for gene, records in self.load_gene_json(self.path,
                                                 file_filter='*.json'):
            print(f"{gene}")
            print(f"# of records = {len(records)}")
            parsed_records = []
            for record in records:
                record = bsoup(str(record), features='lxml')
                parsed_data = ClinVar_Parse_Record(record).get_record_data()
                parsed_records.append(parsed_data)
            gene_records[gene] = parsed_records
        if self.return_data:
            return gene_records
        else:
            # TODO
            pass

    def load_gene_json(self, path, file_filter=None):
        """generator function to load files based on file_filter
        and yields gene, data
        (assumes json encodes a dict with {gene: clinvar_data}

        Arg:
            file_filter = str to search for directory with glob

            path = path to directory with target files
        """
        for filepath in Path(path).glob(file_filter):
            print(f"loading {filepath}")
            if 'parse' in filepath.stem:
                # extracts gene name from "gene_parse.json"
                gene = filepath.stem[:-6]
            else:
                # extracts gene name from "gene.json"
                gene = filepath.stem
            if gene in self.gene_list:
                # checks for gene.json with gene in panel list
                if not self.suppress_output:
                    print(f"loading: {filepath}")
                with open(filepath) as file_:
                    if self.decoder == 'sjson':
                        clinvar_data = self.sjson_decode(
                            file_,
                            decode_type=self.decode_type)
                    else:
                        clinvar_data = self.json_decode(
                            file_,
                            decode_type=self.decode_type)
                yield gene, clinvar_data
            else:
                continue

    @staticmethod
    def read_column(file_path, index=0):
        """simple function to read a line in file as list
        default row 0 (index argument)
        """
        with open(file_path) as f:
            data = []
            reader = csv.reader(f)
            data = [row[index] for row in reader]
        return data
