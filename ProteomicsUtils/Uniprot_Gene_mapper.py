import pprint
import requests
import re, os
import pandas as pd
from pandas.compat import StringIO
#personal modules
from ProteomicsUtils import FileHandling
from ProteomicsUtils.LoggerConfig import logger_config


logger = logger_config(__name__)
logger.info("Import OK")

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def pass_and_retrieve(url, params):
    contact = "dezerae53@hotmail.com" # Please set your email address here to help us debug in case of problems.
    headers={'User-Agent': 'Python %s' % contact}
    request = requests.get(url, params=params, headers=headers, stream=True)

    rawData = request.content
    decoded_data = StringIO(rawData.decode("utf-8"))
    df = pd.read_csv(decoded_data, sep="//t|//n|\\t")
    #logger.info(f'Retrieved dataframe: {df}')
    return df


def protein_info(proteins):
    url = 'https://www.uniprot.org/uniprot/'
    df_list = []
    for protein in proteins:
        params = {'query': protein, 'format': 'tab', 'organism':'mouse', 'reviewed':'yes', 'columns': 'id,entry_name,genes,protein_names,go,comment(FUNCTION)',}
        protein_df = pass_and_retrieve(url, params)
        df_list.append(protein_df)
    return pd.concat(df_list)

def gene_mapper(proteins):
    url = 'http://www.uniprot.org/uploadlists/'

    if len(proteins)<1000:
        proteins = (' ').join(proteins)
        params = {'from':'ACC', 'to':'ID', 'format':'tab', 'query': proteins, 'columns': 'id,entry_name,genes,protein_names,go,comment(FUNCTION)',}
        df = pass_and_retrieve(url, params)
    else:
        dfs = []
        for x in list(chunks(proteins, 1000)):
            protein_chunk = (' ').join(x)
            logger.info(f'protein_chunk: {len(protein_chunk)}')
            params = {'from':'ACC', 'to':'ID', 'format':'tab', 'query': protein_chunk, 'columns': 'id,entry_name,genes,protein_names,go,comment(FUNCTION)',}
            df_chunk = pass_and_retrieve(url, params)
            dfs.append(df_chunk)
            logger.info(f'df_chunk: {df_chunk.shape}')
    return pd.concat(dfs, ignore_index=True)


def main(input, output_path, col_name):

    #Checking input type: expects df or file path
    if isinstance(input, pd.DataFrame):
        logger.info(f"Input detected as DataFrame.")
        dataframe = input
    elif os.path.isfile(input):
        logger.info(f"Collecting proteins from {input}...")
        dataframe = pd.read_excel(input)
    else:
        logger.info(f"Incorrect input format detected. Please pass full file path, or a dataframe as input.")

    proteins = dataframe[col_name].tolist()
    logger.info(f"{len(proteins)} proteins collected from {col_name}")

    logger.info("Fetching geneIDs...")
    gene_mapped = gene_mapper(proteins)
    logger.info(f"Proteins mapped: {gene_mapped}")
    #collect the new identifiers
    proteins_mapped = list(gene_mapped['To'])
    logger.info(f"Gathering gene info...")
    gene_info = protein_info(proteins_mapped)

    logger.info(f"Gene info gathered: {gene_info}")
    FileHandling.df_to_excel(output_path=output_path+'Gene_Mapping_', sheetnames=['Gene ID Mapping', 'Gene Info'], data_frames=[gene_mapped, gene_info])
    logger.info(f"Gene info saved to {output_path}")

    return gene_info


if __name__ == "__main__":
    #default to test file if args not parsed
    input = 'C:/Users/dezer_000/Documents/App_Dev_Projects/ProteomicsUtils/test_data/test_data_Compiled.xlsx'
    output_path = 'C:/Users/dezer_000/Downloads/'
    col_name = 'Accession'
    main(input, output_path, col_name)
