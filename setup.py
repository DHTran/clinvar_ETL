import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name='clinvar_ETL',
    version='0.1.0',
    author='David Tran',
    author_email='davhutra@gmail.com',
    description='Pull and parse ClinVar Variation records',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/DHTran/clinvar_ETL',
    project_urls = {
        "Bug Tracker": "https://github.com/DHTran/clinvar_ETL/issues"
    },
    license='MIT',
    packages=['clinvar_ETL'],
    install_requires=[
        'pandas', 'simplejson', 'beautifulsoup4', 'biopython',
        'lxml', 'google-api-python-client', 'google-auth-oauthlib',
        'google-auth-httplib2', 'python-dotenv'
    ],
)