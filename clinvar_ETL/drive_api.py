"""Create Google sheet from csv file.

TODO: update redirects, does not work currently
"""
import os.path
import pickle
import time
import socket
import sys
from clinvar_ETL.clinvar_utilities import PARSE_CSVS_PATH
from googleapiclient.http import MediaFileUpload
from googleapiclient.errors import HttpError
from googleapiclient.discovery import build
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request
from pathlib import Path
from dotenv import load_dotenv

load_dotenv()
secrets_path = os.getenv('secrets_path')
sys.path.insert(1, secrets_path)
DRIVE_SECRETS_PATH = Path(secrets_path)
CLINVAR_FOLDER = Path('/Volumes/GoogleDrive/My Drive/ClinVar datapulls')

from drive_secrets import SCOPES, CLINVAR_DATAPULLS_FOLDER


class ServiceMixin:
    """Create Google API service.
    """
    def get_service(self, api='sheets', version='v4'):
        """gets tokens using credentials.json if tokens don't already
        exist.  Creates service with build, default api = 'sheets',
        version = 'v4'.
        """
        creds = None
        # The file token.pickle stores the user's access and refresh
        # tokens, and is created automatically when the authorization
        # flow completes for the first time.
        if os.path.exists(DRIVE_SECRETS_PATH/'token.pickle'):
            print("token.pickle exists")
            with open(DRIVE_SECRETS_PATH/'token.pickle', 'rb') as token:
                creds = pickle.load(token)
        # If there are no (valid) credentials, let the user log in.
        if not creds or not creds.valid:
            print("no creds or not valid")
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


class Create_Sheets(ServiceMixin):
    """use Google api to create and update Google sheets using
    ClinVar data.  Takes list of genes, creates filenames with
    gene and current date (day-month-year)

    Args:
        path = path to folder containg csvs
        overwrite (bool) = if True will overwrite sheet files
        folder_id = Google Drive folder to place sheet files
    """

    def __init__(self, path=None, overwrite=False):
        if path:
            self.path = path
        else:
            self.path = PARSE_CSVS_PATH
        self.overwrite = overwrite

    def __repr__(self):
        repr_string = (f"""
        path = {self.path}
        overwrite = {self.overwrite}
        SCOPES = {SCOPES}
        """)
        return repr_string

    def sheets_from_csvs(self):
        """creates sheets from csv files
        """
        # prevent timeout error
        socket.setdefaulttimeout(150)
        service = self.get_service(api='drive', version='v3')
        csv_paths = self.load_filenames(self.path, '*.csv')

        def filter_existing(csv_paths):
            """removes existing sheets from list of paths
            """
            existing_sheets = self.load_filenames(
                CLINVAR_FOLDER, '*.gsheet', return_paths=False)
            # filter out existing sheets
            return [item for item in csv_paths
                    if item[0] not in existing_sheets]

        if not self.overwrite:
            csv_paths = filter_existing(csv_paths)
        for filename, csv_path in csv_paths:
            self.create_sheet(service, filename, csv_path)
            time.sleep(0.1)

    def create_sheet(self, service, filename, csv_path):
        """Create Google sheets from csv files using Drive API.

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


class Update_Sheets(ServiceMixin):
    """Update Google sheets with custom column widths.

    Args:
        api = API service to use, default api is 'sheet'
        version = Verson of the API, for 'sheet' it's 'v4
        date_text = text string to filter on for modification. e.g.
            sheets are named gene_03-27-21, so filter on '03-27-21'
        folder_id = id of Google folder containing sheets
    """

    def __init__(self, date_text, folder_id=None,
                 api='sheet', version='v4'):
        self.api = api
        self.version = version
        self.date_text = date_text
        if folder_id:
            self.folder_id = folder_id
        else:
            self.folder_id = CLINVAR_DATAPULLS_FOLDER

    def __repr__(self):
        repr_string = (f"Update Sheets\
        \n api = {self.api} \
        \n version = {self.version} \
        \n date_text = {self.date_text}")
        return repr_string

    def update_sheets(self, redo=False):
        drive_service = self.get_service(api='drive', version='v3')
        file_ids = self.get_file_ids(drive_service, self.date_text,
                                     self.folder_id)
        sheet_service = self.get_service(api='sheets', version='v4')
        file_sheet_ids = self.get_sheet_ids(sheet_service, file_ids)
        remainder = self.custom_column_widths(sheet_service,
                                              file_sheet_ids)
        # tries 1 more time with items in remainder
        if redo:
            time.sleep(20)
            try:
                remainder = self.custom_column_widths(
                    sheet_service, remainder)
            except HttpError as e:
                print(f"An error occurred: {e}")
        return remainder

    def get_file_ids(self, drive_service, date_text, folder_id):
        """Pull spreadsheetIds of spreadsheets in folder (using folderID)
        returns list of file ids.

        checks if spreadsheet names contains 'date_text'

        e.g. clinVar data named 'ATM_03-26-2021'
        """

        def sort_tuples(tuple_list):
            tuple_list.sort(key = lambda x: x[0])
            return tuple_list

        spreadsheet_ids = []
        page_token = None
        while True:
            try:
                response = drive_service.files().list(
                    q=(f"mimeType='application/vnd.google-apps.spreadsheet' \
                    and parents in '{folder_id}'"),
                    spaces='drive',
                    fields='nextPageToken, files(id, name, modifiedTime)',
                    pageToken=page_token).execute()
                counter = 0
                for file in response.get('files', []):
                    counter += 1
                    name = file.get('name')
                    id_ = file.get('id')
                    modified_date = file.get('modifiedTime')
                    if date_text in file.get('name'):
                        spreadsheet_ids.append((name, id_, modified_date))
                page_token = response.get('nextPageToken', None)
                if page_token is None:
                    print(f"file count = {counter}")
                    break
            except HttpError as e:
                print(f"An error occurred: {e}")
                break
        spreadsheet_ids = sort_tuples(spreadsheet_ids)
        return spreadsheet_ids

    def get_sheet_ids(self, sheet_service, file_ids):
        """Use spreadsheetIds to pull sheet_ids (first sheet) from Google sheets.

        Args: file_ids = list of tuples
            ('name', 'spreadsheet_ids', 'modifiedTime')

        Returns: list of tuples ('name', 'spreadsheetId', 'sheet_id')
        """
        file_sheet_ids = []
        for item in file_ids:
            sheet_data = (
                sheet_service.spreadsheets().
                get(spreadsheetId=item[1]).execute()
                )
            # to prevent going over quota
            time.sleep(5)
            sheets = sheet_data.get('sheets', '')
            # title = sheets[0].get("properties", {}).get("title", "Sheet1")
            sheet_id = sheets[0].get("properties", {}).get("sheetId", 0)
            new_tuple = (item[0], item[1], item[2], sheet_id)
            print(new_tuple)
            file_sheet_ids.append(new_tuple)
        return file_sheet_ids

    # TODO add code to filter out files that have been modified recently
    def custom_column_widths(self, sheet_service, file_sheet_ids):
        """Update specific columns to custom pixel widths.

        Arg: ids_list = list of tuples
            ('name', 'spreadsheetId', 'modifiedTime', 'sheet_id')

        Return: remainder (list of tuples that had HttpError, i.e.
        failed to update)
        """
        remainder = []
        for name, spreadsheet_id, modified_time, sheet_id in file_sheet_ids:
            update_requests = {
                "requests": [
                    {"updateDimensionProperties": {
                        "range": {"sheetId": sheet_id,
                                  "dimension": "COLUMNS",
                                  "startIndex": 1,
                                  "endIndex": 2},
                        "properties": {"pixelSize": 250},
                        "fields": "pixelSize"}, },
                    {"updateDimensionProperties": {
                        "range": {"sheetId": sheet_id,
                                  "dimension": "COLUMNS",
                                  "startIndex": 2,
                                  "endIndex": 3},
                        "properties": {"pixelSize": 100},
                        "fields": "pixelSize"}},
                    {"updateDimensionProperties": {
                        "range": {"sheetId": sheet_id,
                                  "dimension": "COLUMNS",
                                  "startIndex": 3,
                                  "endIndex": 5},
                        "properties": {"pixelSize": 400},
                        "fields": "pixelSize"}},
                    ]
            }
            try:
                request = (
                    sheet_service.spreadsheets().
                    batchUpdate(spreadsheetId=spreadsheet_id,
                                body=update_requests)
                    )
                time.sleep(5)
                response = request.execute()
            except HttpError as e:
                print(f"error: {e}")
                print(f"error on {name}")
                errored = (name, spreadsheet_id, modified_time, sheet_id)
                remainder.append(errored)
        return remainder
