import csv
from pathlib import Path

# constants
DATABASE = 'clinvar'
GLAUCOMA_GENES = 'glaucoma_genes.csv'
BLOCKLIST_PATH = Path('clinvar_ETL/blocklist_pmids.csv')
TEST_GENE = []
MULTIGENES = ['RAD51D']
HOME_PATH = Path.home()
ETL_DATAFILES = Path(HOME_PATH/'datafiles/clinvar_ETL_datafiles')
CLINVAR_DOWNLOADS = ETL_DATAFILES/'ClinVar_downloads'
DATAPULL_JSONS_PATH = ETL_DATAFILES/'datapull_jsons'
GENOME_FILES_PATH = Path(HOME_PATH/'datafiles/genome_datafiles')
PARSE_JSONS_PATH = ETL_DATAFILES/'parse_jsons'
PARSE_CSVS_PATH = ETL_DATAFILES/'parse_csvs'
TEST_RECORDS_PATH = ETL_DATAFILES/'test_records'
TESTS_PATH = ETL_DATAFILES/'tests'
TEST_JSONS_PATH = ETL_DATAFILES/'test_jsons'
TEST_RECORDS_IDS = []
TEST_RECORDS_PATHS_LIST = [f"{id_}_record" for id_ in TEST_RECORDS_IDS]
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
