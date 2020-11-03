from pathlib import Path

# constants
DATABASE = 'clinvar'
CANCER_GENE_SHEET = 'cancer_genes.csv'
CARDIO_GENE_SHEET  = 'cardio_genes.csv'
ACMG59_GENE_SHEET  = 'ACMG59_genes.csv'
# GREM1 only has 12 records
TEST_GENE = ['GREM1']
MULTIGENES = ['RAD51D']
HOME_PATH = Path.home()
DATAFILES_PATH = Path(HOME_PATH/'datafiles/clinvar_ETL_datafiles')
CLINVAR_DOWNLOADS = DATAFILES_PATH/'ClinVar_downloads'
DATAPULL_JSONS_PATH = DATAFILES_PATH/'datapull_jsons'
PARSE_JSONS_PATH = DATAFILES_PATH/'parse_jsons'
PARSE_CSVS_PATH = DATAFILES_PATH/'parse_csvs'
BLOCKLIST_PATH = DATAFILES_PATH/'blocklist_pmids.csv'
TESTS_PATH = DATAFILES_PATH/'tests'
TEST_RECORDS_IDS = [55634, 55628, 127898, 464706, 37003, 228361]