import csv
import simplejson as sjson
import re
import pandas as pd
# patterns module contains regex and text matches
import patterns
from bs4 import BeautifulSoup as bsoup
from datetime import date
from collections import Counter, namedtuple, defaultdict
from pathlib import Path
from clinvar_ETL_constants import DATAFILES_PATH
from clinvar_ETL_constants import DATAPULL_JSONS_PATH, BLOCKLIST_PATH
from clinvar_ETL_constants import PARSE_JSONS_PATH
from clinvar_datapull import ClinVar_Datapull as clinvar_dp

# updated parse class that reads bsoup objects
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
        variation_summary_data = self.variation_summary_data(variant_data)
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
                            pmids.append(item.id.string.strip())
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
                submitter_data['condition'] = disease.strip()
                yield submitter_data
        assertion_results = (
            list(get_submitter_data(clinical_assertions))
            )
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
            submitter['classification'] for submitter in classification_data)
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

                Arg: data = variant_data dict
                """
                classification_submitters = defaultdict(list)
                for submission in submission_data:
                    submitter = submission.get('submittername')
                    classification = submission.get('classification')
                    classification_submitters[classification].append(submitter)
                return classification_submitters

            def clean_classifications(classifications):
                """
                Remove python characters and group

                Args:
                    classifications (list of dicts):
                    [{'submittername': classification}, ...]

                Return:
                    string = 'classifications': submitter1, submitter2...
                """
                result_string = ''
                for classification, submitters in classifications.items():
                    class_submitters = ''
                    classification = classification.string+": "
                    submitters = ClinVar_Parse_Record.clean_string(submitters)
                    class_submitters = "*  " + classification + submitters + ",  "
                    result_string = result_string + class_submitters
                return result_string

            classification_submitters = dict(group_submitters(data['submission_data']))
            cleaned_classifications = (
                clean_classifications(classification_submitters))
            return cleaned_classifications

        def get_submitter_notes(submission_data):
            """method to get 'submitter': submitter notes
            from submission data
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

        summary_dict = {}
        summary_dict['classification_count'] = (
            classification_count(data['submission_data']))
        summary_dict['classifications'] = group_classifications(data)
        summary_dict['comments'] = get_submitter_notes(data['submission_data'])
        keys_to_add = [
            'variation_id', 'name', 'date_updated', 'num_submitters']
        for key in keys_to_add:
            summary_dict[key] = data[key]
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

class Read_jsons:
    """
    class to read variation data pulled from ClinVar and stored as json

    Arg:
        path = path to save parse jsons. Default is DATAPULL_JSONS_PATH

        test_flag(bool):  if test_flag = True, returns record(s) and
        count instead of saving to json using gene.json file name.
        Default = False

        overwrite(bool):  overwrite = True saves over existing file,
            Default = False

        genes(list): user can input custom list of genes, default:  pulls
        list of cancer and cardio genes (using imported ClinVar_Datapull class)

        suppress_out(bool):  True =  prints smaller number of print
            statements that log progress. Default = True
    """
    def __init__(self, path=None, overwrite=False, test_flag=False,
                 genes=None, suppress_output=True):
        if path:
            self.path = path
        else:
            self.path = DATAPULL_JSONS_PATH
        self.overwrite = overwrite
        self.test_flag = test_flag
        self.blocklist_pmids = self.read_column(BLOCKLIST_PATH)
        if genes is None:
            # ClinVar_Datapull by default creates genes list from
            # cancer and cardio lists
            # for future provide option to select genes lists
            genes = clinvar_dp().genes
            self.gene_list = genes
        else:
            if not isinstance(genes, list):
                genes = [genes]
            self.gene_list = genes
        self.suppress_output = suppress_output

    def __repr__(self):
        repr_string = (f"""
        path {self.path},
        test_flag is {self.test_flag}
        overwrite is {self.overwrite}
        gene_list = {self.gene_list}
        blocklist_pmids = {self.blocklist_pmids}
        """)
        return repr_string

    def get_datapull(self):
        """read json file storing variation records by gene
        """
        for gene, clinvar_data in self.load_gene_json(
            self.path, file_filter='*.json'):
            ids_records = clinvar_data[gene]
        return ids_records

    def load_gene_json(self, path, file_filter=None):
        """generator function to load files based on file_filter
        and yields gene, data
        (assumes json encodes a dict with {gene: clinvar_data}

        Arg:
            file_filter = str to search for directory with glob

            path = path to directory with target files
        """
        for file_ in Path(path).glob(file_filter):
            if 'parse' in file_.stem:
                # extracts gene name from "gene_parse.json"
                gene = file_.stem[:-6]
            else:
                # extracts gene name from "gene.json"
                gene = file_.stem
            if gene in self.gene_list:
                # checks for gene.json with gene in panel list
                if not self.suppress_output:
                    print(f"loading: {file_}")
                with open(file_) as f:
                    clinvar_data = sjson.load(f)
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