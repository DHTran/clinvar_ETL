import os.path
import pickle
import time
import socket
import sys
from clinvar_ETL_constants import DATAFILES_PATH
from googleapiclient.http import MediaFileUpload
from googleapiclient.discovery import build
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request
from pathlib import Path
sys.path.insert(1, '/Users/dht/secrets')
from drive_secrets import SCOPES, CLINVAR_DATAPULLS_FOLDER

DRIVE_SECRETS_PATH = Path('/Users/dht/secrets')
# os path to ClinVar Datapulls folder
CLINVAR_FOLDER = Path('/Users/dht/Google Drive/ClinVar datapulls')


class Drive_Api:
    """use Google api to create and update Google sheets using
    ClinVar data.  Takes list of genes, creates filenames with
    gene and current date (day-month-year)

    """

    def __init__(self, path=None, overwrite=False):
        if path:
            self.path = path
        else: 
            self.path = DATAFILES_PATH/'clinvar_datapull_datafiles/parse_csvs'
        self.overwrite = overwrite
        
    def __repr__(self):
        repr_string = (f"""
        path = {self.path}
        overwrite = {self.overwrite}
        SCOPES = {SCOPES}
        folder = {CLINVAR_DATAPULLS_FOLDER}
        """)
        return repr_string

    def get_service(self, api, version):
        """gets tokens using credentials.json if tokens don't already
        exist.  Creates service with build, default api = 'sheets', 
        version = 'v4'.  
        """
        creds = None
        # The file token.pickle stores the user's access and refresh 
        # tokens, and is created automatically when the authorization 
        # flow completes for the first time.
        if os.path.exists(DRIVE_SECRETS_PATH/'token.pickle'):
            print(f"token.pickle exists")
            with open(DRIVE_SECRETS_PATH/'token.pickle', 'rb') as token:
                creds = pickle.load(token)
        # If there are no (valid) credentials, let the user log in.
        if not creds or not creds.valid:
            print(f"no creds or not valid")
            if creds and creds.expired and creds.refresh_token:
                print("refresh token")
                creds.refresh(Request())
            else:
                print("loading from credentials.json")
                flow = InstalledAppFlow.from_client_secrets_file(
                    DRIVE_SECRETS_PATH/'credentials.json', SCOPES)
                creds = flow.run_local_server(port=0)
            # Save the credentials for the next run
            with open(DRIVE_SECRETS_PATH/'token.pickle', 'wb') as token:
                pickle.dump(creds, token)
        service = build(api, version, credentials=creds)
        return service

    def create_sheets_from_csvs(self):
        """creates sheets from csv files
        """
        # prevent timeout error
        socket.setdefaulttimeout(150)
        service = self.get_service(api='drive', version='v3')
        csv_paths = self.load_filenames(self.path, '*.csv')

        def filter_existing(csv_paths):
            existing_sheets = self.load_filenames(
                CLINVAR_FOLDER, '*.gsheet', return_paths=False)
            # filter out existing sheets
            return [
                item for item in csv_paths 
                if item[0] not in existing_sheets
                ]

        if not self.overwrite:
            csv_paths = filter_existing(csv_paths)
        for filename, csv_path in csv_paths:
            self.create_sheet(service, filename, csv_path)
            time.sleep(0.1)
        
    def create_sheet(self, service, filename, csv_path):
        """creates sheets from csv files using Drive API
        Takes source_name = filename of csv (gene_month-day-year.csv)
        Creates sheet with name gene_month-day-year

        """
        sheet_name = filename
        sheet_mimetype = 'application/vnd.google-apps.spreadsheet'
        csv_mimetype = 'text/csv'
        metadata = {
            'name': sheet_name,
            'mimeType': sheet_mimetype,
            'parents': [CLINVAR_DATAPULLS_FOLDER],
        }
        media = MediaFileUpload(csv_path, mimetype=csv_mimetype)
        sheet = service.files().create(
            body=metadata, media_body=media).execute()
        if sheet:
            print(f"""
            Imported {csv_path} to {sheet_name}""")

    def load_filenames(self, path, file_filter=None, return_paths=True):
        """load csv files as list of (filename, csv_path)
        """
        file_paths = []
        for file_ in Path(path).glob(file_filter):
            filename = file_.stem
            file_path = file_
            if return_paths:
                results = (filename, file_path)
            else:
                results = filename
            file_paths.append(results)
        return file_paths