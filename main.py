"""CLI for clinvar_ETL.  Runs clinvar datapull, validate, parse, to_csv,
create_sheets, or update_sheets. Note create_sheets and update sheets are
disabled (10-2022).

Arguments
--datapull:  pull data from ClinVar, saves as gene.json
--validate: compares the list of variation ids from a query to what is pulled
  then offers option to add missing ids via a another datapull
--parse:  parse a ClinVar datapull, saves as gene_parse.json
--to_csv: converts parses to csvs and saves as gene_MM_DD_YYYY (where
  MM_DD_YYYY is the date of the conversion)
--create_sheets: creates Google sheets from csv files
--update_sheets: updates formatting to be more user-readable

--check: prints arguments to confirm, without running
--overwrite (-0): overwrites files during saves
--genes (-g): list of genes to query, e.g. 'APC, BRCA1, NBN'
--test (-t): limit query to 10 files or records
--panel (-p): can select panel of genes, e.g. 'glaucoma'

Example: python main.py --datapull --panel 'glaucoma' -o -t will run a
the datapull.  -o indicates overwrite files, -t limits run to 10
file or records.

"""
import argparse
from Bio import Entrez
from pathlib import Path
from clinvar_datapull import ClinVar_Datapull, Validate_Datapull
from clinvar_parse import Read_ClinVar_Jsons
# from drive_api import Create_Sheets, Update_Sheets
# max number of retries after failed HTTP failures
Entrez.max_tries = 5
eutils_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
# secrets_path = Path("~/secrets")

def main():
    args = None
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        '-o', '--overwrite', help="if supplied overwrite saved files",
        action="store_true")
    parser.add_argument(
        '-g', '--genes', type=lambda x: [gene.strip().upper() for gene in
                                         x.split(',')],
        default=None, help='Gene(s) to query provided, \
        e.g. "APC, MLH1, RET", default=None. If None supplied will use query \
        genes on glaucoma panel. Note query is case-insensitive')
    parser.add_argument(
        '-p', '--panel', type=str, default=None, help="Which gene panel to \
        choose: 'glaucoma'",
        choices=['glaucoma', 'avagen'])
    parser.add_argument(
        '-s', '--save_path', type=str, default=None, help="Path to save file \
            folder")
    parser.add_argument(
        '-t', '--test', action="store_true",
        help="Limit query to 10 files/records")
    parser.add_argument(
        '--datapull', action="store_true", help="Run datapull")
    parser.add_argument(
        '--validate', action="store_true", help="Run validation on datapull")
    parser.add_argument(
        '--parse', action="store_true", help="Run parse of ClinVar data")
    parser.add_argument(
        '--to_csv', action="store_true", help="Convert parse.jsons to CSVs")
    # Drive.api needs updating
    # parser.add_argument(
    #    '--create_sheets', action="store_true",
    #    help="Convert CSV to Google sheet")
    # parser.add_argument(
    #    '--update_sheets', action="store_true", help="Update Sheet formatting")
    parser.add_argument(
        '--date', type=str, default=None, help="Date to filter, MM-DD-YYYY")
    parser.add_argument(
        '--check', action="store_true", help="Check arguments, don't run")
    args = parser.parse_args()
    if args.check:
        print(args)
    elif args.datapull:
        data = ClinVar_Datapull(
            overwrite=args.overwrite, test_flag=args.test,
            gene_panel=args.panel, genes=args.genes, path=args.save_path)
        print(data)
        data.get_records()
    elif args.validate:
        validate = Validate_Datapull(genes=args.genes,
                                     gene_panel=args.panel, save_json=True)
        print(validate)
        validate.check_and_update()
    elif args.parse:
        read_datapulls = Read_ClinVar_Jsons(
            return_data=False, test_flag=args.test,
            overwrite=args.overwrite, genes=args.genes)
        print(read_datapulls)
        read_datapulls.parse_datapull_jsons()
    elif args.to_csv:
        read_datapulls = Read_ClinVar_Jsons(
            return_data=False, overwrite=args.overwrite, genes=args.genes)
        print(read_datapulls)
        read_datapulls.parse_to_csv()
    # elif args.create_sheets:
    #    drive = Create_Sheets()
    #    print(drive)
    #    drive.sheets_from_csvs()
    # elif args.update_sheets and args.date:
    #    update = Update_Sheets(args.date)
    #    print(update)
    #    update.update_sheets()


if __name__ == "__main__":
    main()
