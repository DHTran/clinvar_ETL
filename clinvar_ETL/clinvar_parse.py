"""Load ClinVar variation data stored as jsons and extracts specific data."""
import csv
import json
import simplejson as sjson
import re
import pandas as pd
from bs4 import BeautifulSoup as bsoup
from datetime import date
from collections import Counter, defaultdict
from pathlib import Path
from clinvar_utilities import DATAPULL_JSONS_PATH, BLOCKLIST_PMIDS
from clinvar_utilities import PARSE_JSONS_PATH, PARSE_CSVS_PATH
from clinvar_utilities import REPUTABLE_SUBMITTERS
from clinvar_datapull import ClinVar_Datapull as clinvar_dp


def clean_string(string):
    """Remove python/coding-specific characters for better presentation."""
    # not flake8
    PYTHON_CHARS = re.compile("[\{\}\'\[\]]")

    if isinstance(string, str):
        pass
    else:
        string = str(string)
    return re.sub(PYTHON_CHARS, '', string)


class XMLs_to_Soup_Mixin:
    """Mixin to convert ClinVar xmls to bsoup objects."""

    def xmls_to_bsoup(self, records):
        """Convert XMLs to list of bsoup objects."""
        souped_records = []
        for record in records:
            soup_record = bsoup(str(record), features='lxml')
            souped_records.append(soup_record)
        return souped_records


class Load_Jsons_Mixin:
    """Mixin to load json files (ClinVar datapulls and parses)."""

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
        """Read json and yield gene, data.

        Generator, asssumes json encodes a dict with {gene: clinvar_data}.

        Arg:
            file_filter: str to filter directory with glob
            genes: list of genes to check before attempting to decode
            path: default is DATAPULLS_PATH

        Yields:
            gene and clinvar variation xml
        """
        if file_filter is None:
            file_filter = '*'
        for filepath in Path(path).glob(file_filter):
            if 'parse' in filepath.stem:
                # extracts gene name from "gene_parse.json"
                gene = filepath.stem[:-6]
            else:
                # extracts gene name from "gene.json"
                gene = filepath.stem
            # parse_name = f"{gene}_parse.json"
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


class ClinVar_Parse_Record:
    """Parse a single ClinVar Variation record and return selected
    information.

    Arg:
        record:  a ClinVar Variation record (as Beautiful Soup4 object)
    """
    CLASSIFICATION_NAMES = ['pathogenic', 'benign']
    TRUSTED_NAMES = ["trusted_" + name for name in CLASSIFICATION_NAMES]
    ALL_NAMES = list(zip(CLASSIFICATION_NAMES, TRUSTED_NAMES))

    def __init__(self, record):
        """Take single record as class attribute."""
        self.record = record

    def __repr__(self):
        """Return record header."""
        record_header = self.record.find('variationarchive').attrs
        repr_string = f"record: {str(record_header)}"
        return repr_string

    def get_record_data(self, gene):
        """Parse specific data for a given ClinVar Variation Record.

        variant-level data: ['variationid', 'start', 'stop',
            'referenceallelevcf', 'numberofsubmissions', 'datelastupdated',
            'alternateallelevcf', 'chr', 'assembly','chrom']

        submission data (per submission)
        - submitter_name, classification, inheritance, condition, pmid(s),
        (submitter comment is not parsed currently, method
         get_submitter_comments)
        """
        variant_data = self.get_variant_data(gene)
        variant_data['submission_data'], variant_data['submissions'] = (
            self.get_assertion_data())
        variation_summary_data = self.get_variation_summary_data(variant_data)
        return variation_summary_data

    def get_variant_data(self, gene):
        """Pull variant level data from xml.

        Note need to use gene argument as "gene" since some variation records
        may have an overlapping or neighboring gene(s) in the xml

        e.g. [Gene]SBF2 also returns SNV records with SBF2-AS1

        Args:
            gene: gene name (passed from datapull.json)

        Returns dict of {variation data keys: values}, e.g.
        'variation_id': id, 'gene': gene ...
        """
        results_dict = {}
        allele_data = ['start', 'stop', 'referenceallelevcf',
                       'alternateallelevcf', 'chr', 'assembly',
                       'chrom']
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
        # some variants (e.g. id=29907) report gene region as start/stop
        # in one of the sequencelocation tags, filter for positionvcf to avoid
        if len(sequencelocation_tags) > 1:
            for result in sequencelocation_tags:
                if result.attrs.get('positionvcf'):
                    variant_data = result.attrs
                else:
                    variant_data = self.record.find_all(
                        "sequencelocation", assembly="GRCh37")[-1].attrs
        else:
            # when VCF-stype data is not present
            variant_data = self.record.find_all(
                "sequencelocation", assembly="GRCh37")[0].attrs
        for item in allele_data:
            results_dict[item] = variant_data.get(item)
        results_dict['variation_id'] = variation_id
        results_dict['name'] = name
        results_dict['gene'] = gene
        results_dict['ref'] = results_dict.pop('referenceallelevcf')
        results_dict['alt'] = results_dict.pop('alternateallelevcf')
        results_dict['date_updated'] = date_updated
        results_dict['num_submitters'] = num_submitters
        # rename key 'chr' to 'chrom' to avoid chr() function name
        results_dict['chrom'] = results_dict.pop('chr')
        results_dict['chrom'] = f"chr{results_dict['chrom']}"
        return results_dict

    def get_assertion_data(self):
        """Parse 'clinicalassertion' XML node to get submission data.

        Returns:
            assertion_results: [submitter: assertion data]
        """
        clinical_assertions = self.record.select("clinicalassertion")
        num_submitters = len(clinical_assertions)
        assertion_results = list(self.get_submitter_data(clinical_assertions))
        assert len(assertion_results) == num_submitters, \
            f"""
            mismatch between # of submitter_data dicts
            {len(assertion_results)} and # of assertions
            {len(clinical_assertions)}"""
        return assertion_results, num_submitters

    def get_submitter_data(self, assertions):
        """Iterate over assertions to pull out submitter data."""
        # conditions = []

        def set_classifications_counter(classification_name):
            """Nested function to set 1 or 0 for given classification_name.
            """
            trusted_string = str("trusted_"+classification_name)
            # a few submitters have submitted classifications as
            # "non-pathogenic"
            if "non-pathogenic" in submitter_classification:
                submitter_data[classification_name] = 0
            elif classification_name in submitter_classification:
                submitter_data[classification_name] = 1
            else:
                submitter_data[classification_name] = 0
            if ((classification_name in submitter_classification) and
                    (any(item in submitter_name for item in select_submitters))):
                submitter_data[trusted_string] = 1
            else:
                submitter_data[trusted_string] = 0

        for assertion in assertions:
            submitter_data = {}
            pmids = []
            select_submitters = REPUTABLE_SUBMITTERS
            # scv_number = (
            #    assertion.find("clinvaraccession").
            #    attrs['accession'])
            submitter_data['datelastupdated'] = (
                assertion.attrs['datelastupdated'])
            submitter_name = (
                assertion.find("clinvaraccession").attrs['submittername'])
            submitter_data['submittername'] = submitter_name
            submitter_classification = (
                assertion.find("description").string.strip().lower())
            submitter_data['classification'] = submitter_classification
            set_classifications_counter('pathogenic')
            set_classifications_counter('benign')
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

    def get_variation_summary_data(self, data):
        """Summarize the parsed clinvar data in a user friendly format.

        Args:
            data: variant_data dict (from get_assertion_data() method)

        TODO: 'condition' data is often given as reference to a disease
        notation in another database (e.g. OMIM). To make a readable disease/
        condition will need to query those databases to retrieve the disease
        name.
        """
        summary_dict = {}
        summary_string, classification_counts = (
            self.classification_count(data['submission_data']))
        summary_dict['classification_count'] = summary_string
        for classification_name, counts in classification_counts.items():
            num_of_classifications = counts[1]
            classification_count_str = classification_name + "_count"
            if num_of_classifications:
                summary_dict[classification_count_str] = num_of_classifications
            else:
                summary_dict[classification_count_str] = 0
        variation_url = (
            f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{data['variation_id']}")
        classifications = f"{self.group_submitters_by_classification(data)}"
        summary_dict['classifications'] = classifications
        summary_dict['variation_url'] = variation_url
        # summary_dict['comments'] = get_submitter_comments(
        #    data['submission_data'])
        summary_dict['pmids'], summary_dict['blocked_pmids'] = (
            self.get_pmids(data['submission_data']))
        # summary_dict['condition'] = data['condition']
        keys_to_add = [
            'variation_id', 'name', 'date_updated',
            'num_submitters', 'start', 'stop', 'ref', 'alt',
            'chrom', 'gene']
        for key in keys_to_add:
            summary_dict[key] = data.get(key)
        return summary_dict

    def classification_count(self, classification_data):
        """Get summary and number of classifications per ClinVar variation
        record.

        Arg:
            classification_data: list of dicts

        Returns:
            summary_string: readable string showing counts of each
              classification (e.g. 'ClinVar summary - benign: 1')
            classification_counts: Counter elements for classifications
              defined in ALL_NAMES list of tuples
        """
        classification_counts = {}
        total_classifications = Counter(
            submitter['classification'] for submitter in classification_data)
        for classification_name, trusted_name in self.ALL_NAMES:
            classification_counts[classification_name] = Counter(
                submitter[classification_name] for submitter
                in classification_data)
            classification_counts[trusted_name] = Counter(
                submitter[trusted_name] for submitter
                in classification_data)
        # remove python chars
        cleaned_total_string = clean_string(dict(total_classifications))
        summary_string = f"ClinVar summary - {cleaned_total_string}"
        return (summary_string, classification_counts)

    def group_submitters_by_classification(self, data):
        """Reorganize classification submission data.

        From {'submittername': 'classification'} to
             {'classification': ['submitter1', 'submitter2'...]}

        Arg: data = variant_data dict containing 'submission' data
        """
        classification_submitters = defaultdict(list)
        for submission in data['submission_data']:
            submitter = submission.get('submittername')
            classification = submission.get('classification')
            classification_submitters[classification].append(submitter)
        return self.clean_classifications(classification_submitters)

    def get_pmids(self, submission_data):
        """Get pmids from submissions.

        Return string of pmids and string of blocklist pmids found
        """
        def create_pmid_string(submission_list):
            """Create string of pmids seperated by comma."""
            pmid_string = ''
            for pmid in submission_list:
                pmid = pmid + ", "
                pmid_string += pmid
            pmid_string = pmid_string[:-2]
            return pmid_string

        pmids_list = []
        blocked_pmids_list = []
        for submission in submission_data:
            for pmid in submission['pmids']:
                if (pmid not in pmids_list) and (pmid not in BLOCKLIST_PMIDS):
                    pmids_list.append(pmid)
                elif pmid in BLOCKLIST_PMIDS:
                    blocked_pmids_list.append(pmid)
                else:
                    pass
        if pmids_list and blocked_pmids_list:
            assert pmids_list != blocked_pmids_list, (
                f"{pmids_list} and {blocked_pmids_list} should not be the same")
        pmid_string = create_pmid_string(pmids_list)
        blocked_string = create_pmid_string(blocked_pmids_list)
        if pmid_string and blocked_string:
            assert pmid_string != blocked_string
        return pmid_string, blocked_string

    def clean_classifications(self, classifications):
        """Remove python characters and group by classification.

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
            submitters = clean_string(submitters)
            if result_string:
                class_submitters = (
                    "\n*  " + classification + submitters + ", ")
            elif not result_string:
                class_submitters = ("*  " + classification + submitters + ", ")
            result_string = result_string + class_submitters
        return result_string

    def get_submitter_comments(self, submission_data):
        """Get comment data from submissions.

        Note: currently not used
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


class Read_Jsons(Load_Jsons_Mixin, XMLs_to_Soup_Mixin):
    """Read variation data pulled from ClinVar (from xml using beautiful
    soup xml parser), store as json.

    Methods:
        parse_datapulls: reads ClinVar xml from datapull jsons and
        stores files as jsons with filename pattern of gene_parse.json

    Args:
        genes(list): list of genes to filter for, default: pulls
                list of cancer and cardio genes (using imported
                ClinVar_Datapull class)
        return_data(bool): if True returns record(s) and count instead
            of saving to json. Default = False
        overwrite(bool): if True saves over existing file,
            Default = False
    """

    def __init__(self, test_flag=False, overwrite=False, return_data=False,
                 genes=None, gene_panel=None):
        """Initialize overwrite and return_data flags, genes and gene_panel
        attributes.

        Also defines Paths for data folders from imported constants.

        If no genes are provided, calls clinvar_datapull class (clinvar_dp) to
        automatically create gene
        """
        self.parse_jsons_path = PARSE_JSONS_PATH
        self.datapull_jsons_path = DATAPULL_JSONS_PATH
        self.csvs_path = PARSE_CSVS_PATH
        self.overwrite = overwrite
        self.return_data = return_data
        self.test_flag = test_flag
        if genes is None:
            # ClinVar_Datapull by default creates genes list from default panel
            genes = clinvar_dp(gene_panel=gene_panel).genes
            self.genes = genes
        else:
            if not isinstance(genes, list):
                genes = [genes]
            self.genes = genes

    def __repr__(self):
        """Print class attributes and blocklist_pmids."""
        repr_string = (f"""
        Read_Jsons class
        datapull_jsons_path = {self.datapull_jsons_path},
        parse_jsons_path {self.parse_jsons_path},
        csvs_path {self.csvs_path},
        return_data is {self.return_data},
        overwrite is {self.overwrite},
        genes = {self.genes},
        blocklist_pmids = {BLOCKLIST_PMIDS},
        test_flag = {self.test_flag}
        """)
        return repr_string

    def parse_datapull_jsons(self):
        """Read json files storing variation records by gene from
        DATAPULL_JSONS_PATH (default), hard-coded filter '*.json'.

        Return: Dict of records, with genes as the key

        TODO: filter for just the genes in 'genes' if given
        """
        gene_records = {}
        for gene, records in self.load_gene_json(
                self.datapull_jsons_path, self.genes, file_filter='*.json'):
            print(f"{gene}")
            print(f"# of records = {len(records)}")
            parsed_records = []
            soup_records = self.xmls_to_bsoup(records)
            if self.test_flag:
                soup_records = soup_records[0:10]
            for record in soup_records:
                parsed_data = ClinVar_Parse_Record(record).get_record_data(gene)
                parsed_records.append(parsed_data)
            if not self.return_data:
                if self.test_flag:
                    filename = f"{gene}_parse_test.json"
                else:
                    filename = f"{gene}_parse.json"
                # store in PARSE_JSONS_PATH by default
                file_path = Path(self.parse_jsons_path/filename)
                if file_path.exists() and not self.overwrite:
                    print("file exists")
                    print(f"overwrite is {self.overwrite} and skipping")
                    continue
                print(f"saving to {file_path}")
                with open(file_path, 'w') as file:
                    clinvar_dp().store(parsed_records, file)
                continue
            gene_records[gene] = parsed_records
        return gene_records

    def parse_to_csv(self, today=None):
        """Load parse jsons and converts to dataframe and saves as csv.

        Args:
            parse_path: path to folder with parse jsons
        """
        column_order = [
            'gene',
            'variation_id',
            'variation_url',
            'name',
            'classifications',
            'classification_count',
            'date_updated',
            'pmids',
            'blocked_pmids',
            'start',
            'stop',
            'ref',
            'alt',
            'chrom',
            'num_submitters',
            'pathogenic_count',
            'trusted_pathogenic_count',
            'benign_count',
            'trusted_benign_count',
            ]

        if today is None:
            date_today = date.today()
            today = date_today.strftime("%m-%d-%Y")
        for gene, data in self.load_gene_json(
                self.parse_jsons_path, genes=self.genes,
                file_filter='*_parse.json'):
            gene_date = gene + "_" + today
            filename = f"{gene_date}.csv"
            savefile = Path(self.csvs_path/filename)
            if savefile.is_file():
                if self.overwrite:
                    pass
                else:
                    print(f"file {filename} exists in folder")
                    print(f"overwrite is {self.overwrite}")
                    print("skipping save")
                    continue
            dataframe = pd.DataFrame(data)
            dataframe = dataframe[column_order]
            sorted_df = dataframe.sort_values(by=['start'])
            sorted_df.to_csv(self.csvs_path/filename, encoding='utf-8',
                             index=False)
            print(f"saving to csv {filename}")

    @staticmethod
    def read_column(file_path, index=0):
        """Use csv.reader on a file to be returned as a list, starting
        each row at a given position (default = 0).

        Args:
            file_path: path to file
            index: position to start converting the reader output to list
        """
        with open(file_path) as f:
            data = []
            reader = csv.reader(f)
            data = [row[index] for row in reader]
        return data

