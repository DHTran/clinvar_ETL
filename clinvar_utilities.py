import csv
from pathlib import Path

# constants
DATABASE = 'clinvar'
# assumes /clinvar_ETL/clinvar_ETL/clinvar_utilities.py location
PARENT_FOLDER = Path().resolve()
PACKAGE_DATAFILES = Path(PARENT_FOLDER/'datafiles')
DATAPULL_JSONS_PATH = PACKAGE_DATAFILES/'datapull_jsons'
PARSE_JSONS_PATH = PACKAGE_DATAFILES/'parse_jsons'
PARSE_CSVS_PATH = PACKAGE_DATAFILES/'parse_csvs'
TEST_RECORDS_PATH = PACKAGE_DATAFILES/'test_records'
BLOCKLIST_PATH = PACKAGE_DATAFILES/'blocklist_pmids.csv'
GLAUCOMA_GENES = PACKAGE_DATAFILES/'glaucoma_genes.csv'
REPUTABLE_SUBMITTERS = [
    'Invitae', 'GeneDx', 'Color', 'Ambry',
    'Partners', 'Counsyl', 'Insight', 'ARUP',
    'Eurofins', 'Prevention Genetics', 'Quest',
    'Athena', 'Baylor Genetics', 'Fulgent genetics',
    'Emory genetics laboratory',
    'ClinGen Glaucoma Variant Curation Expert Panel, ClinGen'
    ]

def read_column(file_path, index=0):
    """simple function to read a line in file as list
    default row 0 (index argument)
    """
    with open(file_path) as f:
        data = []
        reader = csv.reader(f)
        data = [row[index] for row in reader]
    return data

BLOCKLIST_PMIDS = read_column(BLOCKLIST_PATH)
