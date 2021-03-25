import csv
from pathlib import Path

# constants
DATABASE = 'clinvar'
CANCER_GENE_SHEET = 'cancer_genes.csv'
CARDIO_GENE_SHEET = 'cardio_genes.csv'
ACMG59_GENE_SHEET = 'ACMG59_genes.csv'
# GREM1 only has 12 records
TEST_GENE = ['GREM1']
MULTIGENES = ['RAD51D']
HOME_PATH = Path.home()
DATAFILES_PATH = Path(HOME_PATH/'datafiles/clinvar_ETL_datafiles')
TEST_RECORDS_PATH = DATAFILES_PATH/'test_records'
CLINVAR_DOWNLOADS = DATAFILES_PATH/'ClinVar_downloads'
DATAPULL_JSONS_PATH = DATAFILES_PATH/'datapull_jsons'
PARSE_JSONS_PATH = DATAFILES_PATH/'parse_jsons'
PARSE_CSVS_PATH = DATAFILES_PATH/'parse_csvs'
BLOCKLIST_PATH = DATAFILES_PATH/'blocklist_pmids.csv'
TESTS_PATH = DATAFILES_PATH/'tests'
TEST_JSONS_PATH = DATAFILES_PATH/'test_jsons'
TEST_RECORDS_IDS = [
    236589, 55634, 55628, 127898, 464706, 37003, 228361, 37734,
    233913, 9342
    ]
TEST_RECORDS_PATHS_LIST = [f"{id_}_record" for id_ in TEST_RECORDS_IDS]

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