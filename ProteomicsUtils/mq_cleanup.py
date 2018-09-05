
import pandas as pd
import numpy as np
import os, re

from ProteomicsUtils import FileHandling, StatUtils, PlotUtils, DataWrangling
from ProteomicsUtils.LoggerConfig import logger_config

logger = logger_config(__name__)
logger.info('Import OK')

pd.set_option('display.max_columns', 500)


def multifile_cleaner(input_folder, output_path, sample_names=None, proteins_file='proteinGroups.txt', peptides_file='peptides.txt', prot_cols=[], pep_cols=[]):

    peptides = pd.read_table(input_folder+peptides_file, sep='\t')
    logger.info(f'Imported peptides from {peptides_file}. {peptides.shape[0]} entries found.')

    cleaned_dfs = {}
    standard_cols = ['Sequence', 'Proteins', 'Gene names', 'Protein names', 'Unique (Groups)', 'Unique (Proteins)', ] + pep_cols

    if sample_names is None:
        logger.info(f'Sample names not set. Collecting all samples.')
        sample_names = [x.replace('Identification type ', '').split('_')[:-1] for x in peptides.columns.tolist() if 'Identification type ' in x]
        sample_names = [('_').join(x) for x in sample_names]
        logger.info(f'Samples detected: {sample_names}')


    for sample in sample_names:
        sample_cols = [col for col in peptides.columns.tolist() if 'Ratio H/L '+sample in col]
        cleaned_dfs[sample] = peptides[standard_cols + sample_cols]
    logger.info(f'Successfully cleaned peptide dataframe.')

    logger.info(f'Collecting proteins')
    proteins = pd.read_table(input_folder+proteins_file, sep='\t')
    logger.info(f'Imported proteins from {proteins_file}. {proteins.shape[0]} entries found.')

    # remove contaminant and reverse proteins
    proteins = proteins[(proteins['Reverse'] != '+') & (proteins['Potential contaminant'] != '+')]
    logger.info(f'Removed contaminant and reverse proteins: {proteins.shape[0]} entries remain.')

    cleaned_prot_dfs = {}
    standard_cols = ['Protein IDs', 'Gene names', 'Protein names', 'Number of proteins', ] + prot_cols

    for sample in sample_names:
        sample_cols = [col for col in proteins.columns.tolist() if sample in col]
        #collect columns of interest
        sample_vals = proteins[standard_cols + sample_cols]
        #collect only proteins with at least one peptide identified in that sample
        sample_reps = sample_vals[[col for col in proteins.columns.tolist() if 'Peptides '+sample in col]].sum(axis=1)
        # collect only proteins which are master proteins
        master_proteins = sample_vals[sample_vals['Number of proteins'] == 1]
        cleaned_prot_dfs[sample] = master_proteins[sample_reps > 0]

    logger.info(f'Successfully cleaned proteins dataframe.')

    logger.info(f'Sorting cleaned data per sample...')
    ## Collecting specific results for each set of samples for further processing
    for sample in sample_names:
        #collect peptide dataframe, rename relevant columns
        pep_dataframe = cleaned_dfs[sample]
        MQ_cols = ['Protein IDs', 'Proteins', 'Protein names'] + [col for col in proteins.columns.tolist() if 'Ratio H/L '+sample in col]
        new_cols = ['ProteinID', 'ProteinID', 'Description'] + [f'Abundance Ratio: ({sample}_{x})' for x in range(1, len(MQ_cols)+1)]
        pep_dataframe.rename(columns=dict(zip(MQ_cols, new_cols)), inplace=True)

        # collect protein dataframe, rename relevant columns
        prot_dataframe = cleaned_prot_dfs[sample]
        prot_dataframe.rename(columns=dict(zip(MQ_cols, new_cols)), inplace=True)

        #save to individual excel spreadsheets
        FileHandling.df_to_excel(output_path+sample+'_MQ_Proteins.xlsx', sheetnames=['Proteins'], data_frames = [prot_dataframe])
        FileHandling.df_to_excel(output_path+sample+'_MQ_Peptides.xlsx', sheetnames=['Peptides'], data_frames = [pep_dataframe])
        FileHandling.df_to_excel(output_path+sample+'_Compiled.xlsx', sheetnames=['Proteins', 'Peptides'], data_frames = [prot_dataframe, pep_dataframe])
        #logger.debug(pep_dataframe.columns.tolist())

    logger.info(f'Proteins and peptides successfully cleaned. Dataframes save to {output_path}.')

    return cleaned_dfs

if __name__ == '__main__':
    input_folder = 'C:/Users/dezer_000/Documents/App_Dev_Projects/ProteomicsUtils/test_data/'

    output_path = 'C:/Users/dezer_000/Documents/App_Dev_Projects/ProteomicsUtils/test_data/'

    multifile_cleaner(input_folder, output_path)
