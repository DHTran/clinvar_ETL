
file_ids = []
page_token = None
while True:
    try:
        response = service.files().list(
            q=(f"mimeType='application/vnd.google-apps.spreadsheet' \
            and parents in '{update.folder_id}'"),
            spaces='drive',
            fields='nextPageToken, files(id, name)',
            pageToken=page_token).execute()
        for file in response.get('files', []):
            # Process change
            print(f"Found file: {file.get('id')}")
            print(f"name: {file.get('name')}")
            file_ids.append(file.get('id'))
        page_token = response.get('nextPageToken', None)
        if page_token is None:
            break
    except HttpError as e:
        print(f"An error occurred: {e}")
        break

- - -


spreadsheetName = "###"  # Please set Spreadsheet name.
sheetName = "###"  # Please set sheet name.

ss = client.open(spreadsheetName)
sheetId = ss.worksheet(sheetName)._properties['sheetId']
body = {
    "requests": [
        {
            "autoResizeDimensions": {
                "dimensions": {
                    "sheetId": sheetId,
                    "dimension": "COLUMNS",
                    "startIndex": 0,  # Please set the column index.
                    "endIndex": 2  # Please set the column index.
                }
            }
        }
    ]
}
res = ss.batch_update(body)


GET https://www.googleapis.com/drive/v3/files?q=mimeType%20!%3D%20%22application%2Fvnd.google-apps.folder%22&fields=files%2Fid&key=[YOUR_API_KEY] HTTP/1.1

GET https://www.googleapis.com/drive/v2/files/folderId/children


Authorization: Bearer [YOUR_ACCESS_TOKEN]
Accept: application/json


results = service.files().list(
    pageSize=10, fields="nextPageToken, files(id, name)").execute()
items = results.get('files', [])