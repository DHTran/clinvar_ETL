import csv
import simplejson as sjson
import re
import pandas as pd
# patterns module contains regex and text matches
import patterns
from datetime import date
from collections import Counter, namedtuple
from pathlib import Path
from clinvar_ETL_constants import GENE_FILES, DATAFILES_PATH
from clinvar_ETL_constants import DATAPULL_JSONS_PATH, BLOCKLIST_PATH
from clinvar_ETL_constants import PARSE_JSONS_PATH
from clinvar_datapull import ClinVar_Datapull as clinvar_dp


class ClinVar_Parse:
    """Parses clinvar variation pages extracting these data:
    allele data: variation id, chrom, start, stop, ref, alt
    submitters, classification, counts of classification types,
    pmids associated with submitters, date updated
    
    ClinVar Variation records inputted as dict of lists:
    {'gene1': [(id, record), (id, record)..], 'gene2': [...]}
    
    default datapull_path folder contains gene.json files

    if test_flag = True, returns record instead of saving to json
    using gene.json file name
    if overwrite = True saves over existing file

    gene_list is used to only load jsons with a filename containing a 
    gene created from ClinVar_Datapull gene list (e.g. cardio genes)

    blocklist_pmids contains list of pmids that are not case or functional
    studies (manually curated)
    """

    def __init__(self, path=None, overwrite=False, test_flag=False):
        if path:
            self.path = path
        else:
            self.path = DATAPULL_JSONS_PATH
        self.overwrite = overwrite
        self.test_flag = test_flag
        self.blocklist_pmids = self.read_column(BLOCKLIST_PATH)
        self.gene_list = (
            clinvar_dp.create_gene_list(GENE_FILES, DATAFILES_PATH))
    
    def __repr__(self):
        repr_string = (f"""
        path {self.path}, 
        test_flag is {self.test_flag}
        overwrite is {self.overwrite}
        gene_list = {self.gene_list}
        blocklist_pmids = {self.blocklist_pmids}
        """)
        return repr_string
        
    def get_variation_data(self):
        """main caller to parse clinvar data from jsons.  
        Loads jsons, calls parsing functions, stores as gene_parse.json
        """
        count_total = []
        # self.load_gene_json yields a json
        for gene, clinvar_data in self.load_gene_json(
                self.path, file_filter='*.json'):
            print(f"parsing records from {gene}")
            ids_records = clinvar_data[gene]
            parse_data, parsed_count = (
                self.parse_clinvar_records(ids_records)
                )
            parse_filename = f"{gene}_parse.json" 
            parse_exists = (
                Path(PARSE_JSONS_PATH/parse_filename).is_file()
                )
            if self.test_flag:
                return parse_data, parsed_count
            elif parse_exists and not self.overwrite:
                print(f"skipping {parse_filename}") 
                continue
            else:
                clinvar_dp.store_data_as_json(
                    parse_data, PARSE_JSONS_PATH, parse_filename)
                gene_parsed_count = (gene, parsed_count)
                count_total.append(gene_parsed_count)
        count_filename = "gene_parsed_count.json"
        # saves json with counts of records per gene
        clinvar_dp.store_data_as_json(
            count_total, PARSE_JSONS_PATH, count_filename)

    def parse_clinvar_records(self, ids_records):
        """loops through Clinvar Variation records (as list of text), 
        pulls specific data items via line position (target index) 
        and regex extraction from each line
        
        ids_records input is a list of dicts [{id:record},...]
        Returns parse_data (list of dicts, each dict = one variation id)

        and count of records per variation id
        """
        print("running parse_clinvar_records")
        parse_data = []
        num_records = len(ids_records)
        Records_Total = namedtuple(
            'Records_Total', ['parsed_records', 'total'])
        for id_record in ids_records:
            variation_id, variation_name, allele_data = None, None, None
            submitters, pmids = [], []
            clinvar_data = {}
            [[id_, record]] = id_record.items()
            id_ = int(id_)
            print(f"record id: {id_}")
            # checks for 'VariationType="copy number', if found skip record
            cnv_present = self.filter_copy_number(record)
            if cnv_present:
                print(f"found cnv, skip")
                continue
            target_indices = self.get_target_indices(record)
            variation_id, variation_name, date_updated = (
                self.get_variation_name_date(target_indices, id_, record)
                )
            allele_data = self.get_allele_data(target_indices, record)
            print(f"allele data in parse_clinvar_records: {allele_data}")
            submitters, classification_count, comments = (
                self.get_submitters_and_submitted_data(
                    target_indices, record)
                )
            pmids = self.get_associated_pmids(target_indices, record)
            clinvar_data = {
                'variation_id': variation_id,
                'date_updated': date_updated,
                'chrom': allele_data.chrom,
                'variant': variation_name,
                'ref': allele_data.ref_seq,
                'alt': allele_data.alt_seq,
                'start': allele_data.start,
                'stop': allele_data.stop,
                'submitters': submitters,
                'classification_count': classification_count,
                'pmids': pmids,
                'number_submitters': len(submitters),
                'number_pmids': len(pmids),
                'comments': len(comments),
            }
            clinvar_data.update(comments)
            parse_data.append(clinvar_data)
        parsed_count = Records_Total(
            parsed_records=len(parse_data), total=num_records)
        print(f"Parsed count {parsed_count}")
        return parse_data, parsed_count
    
    def load_gene_json(self, path, file_filter=None):
        """loads files bases on file_filter and yields gene, data
        (assumes json encodes a dict with {gene: clinvar_data}
        """
        print("  in load_gene_json")
        print("  filtering in {path} for {file_filter}")
        for file_ in Path(path).glob(file_filter):
            if 'parse' in file_.stem:
                # from parse json
                gene = file_.stem[:-6]
            else:
                # from datapull json
                gene = file_.stem
            if gene in self.gene_list:  
                # checks for gene.json with gene in panel list
                print(f"loading: {file_}")
                with open(file_) as f:
                    clinvar_data = sjson.load(f)
                yield gene, clinvar_data
            else: 
                continue
                
    def filter_copy_number(self, record):
        """filters out CNVs by find of 'copy number'
        """
        assert isinstance(record, list)
        test_string = str(record)
        cnv = test_string.find('VariationType="copy number')
        if cnv > 0:
            return True
        else:
            return False

    def get_target_indices(self, record):
        """Finds the indices of target text in Clinvar record 
        (originally elements in variation page xml, as list of strings)
        
        Each target text is item in patterns.INDEX_KEYS. Record is 
        tuple of ncbi_id, record text (id, record)
        """
        print("getting target indices")
        target_indices = {key: None for key in patterns.INDEX_KEYS}
        for key in target_indices.keys():
            # VCF_CONDITIONS is tuple with two potential matches 
            # when no VCF data is present, must match both targets in line
            if isinstance(key, tuple):
                target_indices[key] = [
                    index for index, line in enumerate(record) 
                    if (key[0] in line) and (key[1] in line)
                ]
            else:
                target_indices[key] = [
                    index for index, line in enumerate(record) 
                    if key in line
                    ]
        return target_indices
    
    def get_variation_name_date(self, target_indices, id_, record):
        """pulls variation id and variant name 
        """ 
        print(" running get_variation_name_date") 
        index = target_indices.get(patterns.VARIATION_ID_TARGET)[0]
        variation_id = int(
            re.match(patterns.VARIATION_ID_MATCH, record[index]).group(1)
            )
        assert id_ == variation_id, "variation_id does not match id_"
        variation_name = (re.match(
            patterns.VARIATION_NAME_MATCH, record[index]).group(1)
            )
        variation_name = (variation_name.replace('&gt;', '>'))
        date_updated = (
            re.match(patterns.DATE_UPDATED_MATCH, record[index]).group(1)
            )
        return variation_id, variation_name, date_updated
    
    def get_allele_data(self, target_indices, record):
        """record is list of strings that compose the clinvar xml.
        target_indices dict provides the indices for each target text
        to pull chrom, start, stop, alt and ref data. If no VCF data, 
        assume no variant mapping
        
        Return as Allele_data namedtuple 
        ('chrom', 'start', 'stop', 'ref_seq', 'alt_seq')
        
        """
        print(" running get_allele_data") 
        Allele_data = namedtuple(
                'Allele_data', 
                ['chrom', 'start', 'stop', 'alt_seq', 'ref_seq'], 
                defaults=(None, None, None, None, None))
        allele_data = None
        index = target_indices.get(patterns.VCF_CONDITIONS)
        if index:
            index = index[0]
            try:
                start = re.match(
                    patterns.START_MATCH, record[index]).group(1)
                stop = re.match(
                    patterns.STOP_MATCH, record[index]).group(1)
            except AttributeError:
                start = ""
                stop = ""
            try: 
                alt_seq = re.match(
                    patterns.ALT_SEQ_MATCH, record[index]).group(1)
                ref_seq = re.match(
                    patterns.REF_SEQ_MATCH, record[index]).group(1)
            except AttributeError:
                # in some cases no position data is given
                alt_seq = ""
                ref_seq = ""
        elif not index:
            # patterns.VCF_CONDITIONS holds two target texts for position
            # if both present use second: 'referenceAlleleVCF', 
            # else use first: 'SequenceLocation Assembly="GRCh37"'
            if len(target_indices.get(patterns.VCF_CONDITIONS[0])) > 1:
                index = target_indices.get(patterns.VCF_CONDITIONS[0])[1]
            else:
                index = target_indices.get(patterns.VCF_CONDITIONS[0])[0]
            start, stop, alt_seq, ref_seq = None, None, None, None
            print(f"no vcf:  {start}, {stop}")
        chrom = re.match(patterns.CHROM_MATCH, record[index]).group(1)
        allele_data = Allele_data(chrom=chrom, start=start, stop=stop, 
                                  ref_seq=ref_seq, alt_seq=alt_seq)
        return allele_data
    
    def get_submitters_and_submitted_data(self, target_indices, record):
        """pulls submitter names and associated data from ClinVar record
        """
        print(f" running get_submitters_and_submitted_data")
        submitters = []
        comments = {}
        classification_count_dict = {
            'pathogenic': 0, 
            'benign': 0, 
            'unknown': 0, 
            'uncertain': 0,
            'risk': 0,
            'not provided': 0,
            'pathologic': 0,
            'untested': 0,
            'drug response': 0,
            'cancer': 0,
            'association': 0,
            'likely pathogenic': 0,
        }
        indices = target_indices.get(patterns.SUBMITTER_TARGET)   

        def get_classification_and_comments(index, record):
            """parses classifications and submitters comments from ClinVar
            records.  Counts number of different classifications in
            records from a gene.  

            the classification and comments positions may vary after 
            the line containing 'Description>' and 'Comment>' respectively

            Finds line index of the classification and comments, then
            matches the classification and comments contents based on 
            regex patterns
            """
            def get_index(partial_record, target_text):
                """iterates through lines of the partial record
                to find the relative position of a target text (whose
                position may vary) then calculates the index in the 
                whole record
                """
                relative_index = next(
                    (index for index, item in enumerate(partial_record) 
                     if target_text in item), None)
                if relative_index:
                    return relative_index + index

            print("  running get_classification_comments") 
            classification_results, comment_results = None, None
            classification_index, comment_index = None, None
            target_texts = ['Description>', 'Comment>']
            # finds indices of the classification and comments
            for target_text in target_texts:
                absolute_index = get_index(record[index:], target_text)
                if absolute_index:
                    if target_text == 'Description>':
                        classification_index = absolute_index
                    elif target_text == 'Comment>':
                        comment_index = absolute_index  
            # gets classification and comments. If no index then 
            # results is default None
            print(classification_index, comment_index)
            if classification_index:
                classification_results = re.match(
                    patterns.BETWEEN_ANGLE_MATCH, record[classification_index]
                    ).group(1)
            if comment_index:
                comment_results = re.match(
                    patterns.BETWEEN_ANGLE_MATCH, record[comment_index]
                    ).group(1)
            return classification_results, comment_results
        
        def update_count(classification, classification_count_dict):
            """takes classification and classification_count_dict
            and updates classification count by one
            """
            for key in classification_count_dict.keys():
                if key in classification.lower():
                    classification_count_dict[key] += 1
                    return classification_count_dict
    
        for index in indices:
            submitter_dict = {}
            classification, submitter_name = None, None
            submitter_name = re.match(
                patterns.SUBMITTER_NAME_MATCH, record[index]).group(1)
            # filter out OMIM submissions
            if submitter_name == 'OMIM':
                continue         
            classification, comments[submitter_name] = (
                get_classification_and_comments(index, record)
                )
            submitter_dict[submitter_name] = classification
            if classification not in classification_count_dict:
                new_dict = {classification: 0}
                classification_count_dict.update(new_dict)
            try:
                classification_count_dict = update_count(
                    classification, classification_count_dict)
            except TypeError:
                print(classification)
                print(classification_count_dict)
            # classification_count_dict[classification] += count
            submitters.append(submitter_dict)
        return submitters, classification_count_dict, comments  
    
    def get_associated_pmids(self, target_indices, record):
        """pulls out pmids on the front page of the variation
        page
        """
        print(" running get_associated_pmids")
        pmids = []
        start = target_indices.get(patterns.PMID_TARGET)[0]
        try:
            end = target_indices.get(patterns.PMID_END)[0]
        except IndexError:
            # indexError if no pmid_end target
            if 'NumberOfSubmissions="0"' in str(record):
                # 0 submissions means no pmids
                return []              
        for index, line in enumerate(record[start:end]):
            try: 
                pmid = re.match(patterns.PMID_MATCH, line).group(1)
                new_pmid = pmid not in pmids
                not_blocklist = pmid not in self.blocklist_pmids
                if new_pmid and not_blocklist:
                    pmids.append(pmid)
            except AttributeError:
                pass
        return pmids
    
    def convert_parse_to_csv(self, today=None):
        """loads parse jsons and converts to dataframe
        and saves as csv
        """
        if today is None:
            date_today = date.today()
            today = date_today.strftime("%m-%d-%Y")
        jsons_path = DATAFILES_PATH/'clinvar_datapull_datafiles/clinvar_parses'
        csvs_path = DATAFILES_PATH/'clinvar_datapull_datafiles/parse_csvs'
        for gene, data in self.load_gene_json(
                jsons_path, file_filter='*_parse.json'):
            dataframe = pd.DataFrame(data)
            sorted_df = dataframe.sort_values(by=['start'])
            gene_date = gene+"_"+today
            filename = f"{gene_date}.csv"
            file_ = Path(csvs_path/filename)
            if file_.is_file():
                if self.overwrite:
                    pass
                else:
                    print(f"file {filename} exists in folder")
                    print(f"overwrite is {self.overwrite}")
                    print("skipping save")   
                    continue   
            sorted_df.to_csv(csvs_path/filename, encoding='utf-8', index=False)
            print(f"saving to csv {filename}")            
        
    def count_pmids(self):
        """takes list of dicts of Clinvar data where 
        [{'variation_id1': id1, 'pmid': [pmid1, pmid2],...},
         {'variation_id2': id2, 'pmid': []}]
        and counts occurrences of a given pmid number
        
        """
        gene_pmids = {}
        jsons_path = (
            DATAFILES_PATH/'clinvar_datapull_datafiles/parse_jsons')
        for gene, parse in self.load_gene_json(
            jsons_path, file_filter='*_parse.json'):
            pmids = []
            for variation in parse:
                [pmids.append(x) for x in variation['pmids']]
            pmids_count = Counter(pmids).most_common()
            gene_pmids[gene] = pmids_count
        return gene_pmids

    @staticmethod
    def read_column(file_path, index=0):
        """simple function to read column of file as list
        default column 0 (index argument)
        """
        with open(file_path) as f:
            data = []
            reader = csv.reader(f)
            data = [row[index] for row in reader]
        return data

