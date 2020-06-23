import pickle
import os.path
from googleapiclient.discovery import build
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request


SCOPES = [
    'https://www.googleapis.com/auth/spreadsheets',
    'https://www.googleapis.com/auth/drive'
]

TEST_SHEET = '1ChnkRb4qrp-YTcpxFGG3VAqJ6sQkxG1PowQapi0UQME'
CLINVAR_DATAPULLS_FOLDER = '1FXCc2DWmzUAJAZjfPNL-K8Jt6cTmYvMW'


def get_service(api='sheets', version='v4'):
    """gets tokens using credentials.json if tokens don't already
    exist.  Creates service with build, default api = 'sheets', 
    version = 'v4'.  
    """
    creds = None
    # The file token.pickle stores the user's access and refresh tokens, and is
    # created automatically when the authorization flow completes for the first
    # time.
    if os.path.exists('token.pickle'):
        with open('token.pickle', 'rb') as token:
            creds = pickle.load(token)
    # If there are no (valid) credentials available, let the user log in.
    if not creds or not creds.valid:
        if creds and creds.expired and creds.refresh_token:
            creds.refresh(Request())
        else:
            flow = InstalledAppFlow.from_client_secrets_file(
                'credentials.json', SCOPES)
            creds = flow.run_local_server(port=0)
        # Save the credentials for the next run
        with open('token.pickle', 'wb') as token:
            pickle.dump(creds, token)
    service = build(api, version, credentials=creds)
    return service

def create_sheet(service, sheet_name):
    """function to create google sheet
    """
    spreadsheet_data = {'properties': 
        {'name': sheet_name, 
        }
    }
    sheet = (
        service.spreadsheets().create(body=spreadsheet_data).excecute()
    )

def create_sheet(service, sheet_name):
    """creates sheets using Drive API
    """
    DRIVE = build.('drive', 'v3', http=creds.authorize(Http()))
    data = {
        'name': sheet_name,
        'mimeType': 'application/vnd.google-apps.spreadsheet',
        'parents' : [CLINVAR_DATAPULLS_FOLDER],
    }
    sheet = DRIVE.files().create(body=data).execute()

def get_google_sheet(id, service):
    """ 
    Retrieve sheet data using OAuth credentials and Google Python API. 
    """
    gsheet = {}
    print (id)
    gsheet['id'] = id
    # Call the Sheets API
    for key, value in cell_data_dict.items():
        sheet = service.spreadsheets()
        time.sleep(1)
        print (sheet)
        sheet_lookup = (
            sheet.values()
            .get(spreadsheetId=id, range=value)
            .execute()
            )
        cell_values = sheet_lookup.get('values')
        if cell_values is not None:
            gsheet[key] = cell_values[0][0]
        elif cell_values is None:
            gsheet[key] = ""
    return gsheet