from pathlib import Path

# constants
DATABASE = 'clinvar'
GENE_FILES = ['cancer_genes.csv', 'cardio_genes.csv']
# GREM1 only has 12 records
TEST_GENE = ['GREM1']
HOME_PATH = Path.home()
DATAFILES_PATH = Path(HOME_PATH/'datafiles')
CLINVAR_DOWNLOADS = DATAFILES_PATH/'ClinVar_downloads'
MULTIGENES = ['RAD51D']
