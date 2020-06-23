# match digits after start= or stop="
START_MATCH = re.compile('^.*start="(\d*)')
STOP_MATCH = re.compile('^.*stop="(\d*)')
# match digits after Chr="
CHROM_MATCH = re.compile('^.*Chr="(\d*)')
# non-greedy match letters between alternateAlleleVCF=" and " 
ALT_SEQ_MATCH = re.compile('^.*alternateAlleleVCF="(\D.*?)"')
# non-greedy match letters between alternateAlleleVCF=" and "
REF_SEQ_MATCH = re.compile('^.*referenceAlleleVCF="(\D.*?)"')
# match digits after VariationID="
VARIATION_ID_MATCH = re.compile('^.*VariationID="(\d*)')
# match any characters in between quotes after target text
VARIATION_NAME_MATCH = re.compile('^.*VariationName="(.*?)"')
DATE_UPDATED_MATCH = re.compile('^.*DateLastUpdated="([\d-]+)')
SUBMITTER_NAME_MATCH = re.compile('^.*SubmitterName="(.*?)"')
# match any characters in between > <
BETWEEN_ANGLE_MATCH = re.compile('^.*>(.*?)<')
# match digits after <ID Source="PubMed"
PMID_MATCH = re.compile('^.*<ID Source="PubMed">(\d*)')

# target texts in variation record to extract data
VARIATION_ID_TARGET = '<ClinVarResult-Set>'
VCF_CONDITIONS = (
    'SequenceLocation Assembly="GRCh37"', 'referenceAlleleVCF')
GRCH37_TARGET = 'SequenceLocation Assembly="GRCh37"'
SUBMITTER_TARGET = '<ClinVarAccession'
# targets list of pmids on front page
PMID_TARGET = 'Type="Clinical significance">'
PMID_END = '<ConditionList>'
INDEX_KEYS = [variation_id_target, vcf_conditions, grch37_target,
                submitter_target, pmid_target, pmid_end]