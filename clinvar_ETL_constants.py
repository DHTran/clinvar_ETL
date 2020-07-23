from pathlib import Path

# constants
DATABASE = 'clinvar'
CANCER_GENES = 'cancer_genes.csv'
CARDIO_GENES = 'cardio_genes.csv'
# GREM1 only has 12 records
TEST_GENE = ['GREM1']
MULTIGENES = ['RAD51D']
HOME_PATH = Path.home()
DATAFILES_PATH = Path(HOME_PATH/'datafiles')
CLINVAR_DOWNLOADS = DATAFILES_PATH/'ClinVar_downloads'
DATAPULL_JSONS_PATH = DATAFILES_PATH/'clinvar_ETL_datafiles/datapull_jsons'
PARSE_JSONS_PATH = DATAFILES_PATH/'clinvar_ETL_datafiles/parse_jsons'
PARSE_CSVS_PATH = DATAFILES_PATH/'clinvar_ETL_datafiles/parse_csvs'
BLOCKLIST_PATH = DATAFILES_PATH/'clinvar_ETL_datafiles/blocklist_pmids.csv'
TESTS_PATH = DATAFILES_PATH/'clinvar_ETL_datafiles/tests'
TEST_RECORDS_IDS = [55634, 55628, 127898, 464706, 37003, 228361]