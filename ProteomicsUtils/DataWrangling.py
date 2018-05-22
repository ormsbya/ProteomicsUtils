"""
General collection of functions for manipulating dataframes, generally to isolate proteins or peptides that fit the criteria of interest.
"""
import numpy as np
import pandas as pd
from scipy import stats
import os
import logging
from ProteomicsUtils.LoggerConfig import logger_config

logger = logger_config(__name__)
logger.info('Import ok')

def quantified_data(data_frame):
    """Collects only the data for which quantification was completed
    Parameters
    ----------
    data_frame : DataFrame
        Contains raw input from proteomics results preprocessed with ProteomeDiscoverer and exported to excel.

    Returns
    -------
    quant_data: DataFrame
        Contains the filtered dataframe, with Average Abundance ratio column appended.
    col_list: list
        list of column headers containing the abundance ratio data
    """
    #collects only data for which quantifiable info was gathered from PD
    quant_data = data_frame[(data_frame['Quan Info'] == 'Unique')]
    #add a new column which is the average Abundance ratio for multiple replicates
    col_list = [col for col in quant_data.columns if 'Abundance Ratio: (' in col]
    logger.debug('Replicate Columns: {}'.format(col_list))
    AvAbun = []
    for index, row in quant_data.iterrows():
        abundance_ratios = row[col_list]
        av_abun = abundance_ratios.mean()
        AvAbun.append(av_abun)
    quant_data['Abundance Ratio (Average)'] = AvAbun
    logger.debug('Quant_data: {}'.format(quant_data.shape))

    return quant_data, col_list



def Unique_Cys_sorter(quant_data):
    """ Collects proteins for which two unique peptides were found, at least one of which contains a cysteine residue.

    Parameters
    ----------
    quant_data : DataFrame
        Contains raw data from preprocessing of proteomics results, generally one includes quantified proteins

    Returns
    -------
    two_unique_cys, cys_pep, non_cys_pep : DataFrames
        seperate dataframes to collect (1) all proteins that fit the criteria, and individually those that contain (2) or do not contain (3) a cysteine residue.

    """
    # create empty dataframes to store filtered data
    two_unique = pd.DataFrame()
    two_unique_cys = pd.DataFrame()
    cys_pep = pd.DataFrame()
    non_cys_pep = pd.DataFrame()

    for x in quant_data['Master Protein Accessions']:
        protein = quant_data.loc[(quant_data['Master Protein Accessions'] == x)]
        if len(protein)>1:
            cys = protein.loc[(protein['Annotated Sequence'].str.contains('C') == True)]
            non_cys = protein.loc[(protein['Annotated Sequence'].str.contains('C') == False)]
            two_unique = two_unique.append(protein, ignore_index = True)

            if len(cys)>=1 and len(non_cys)>=1:
                two_unique_cys = two_unique_cys.append(protein, ignore_index = True)
                cys_pep = cys_pep.append(cys, ignore_index = True)
                non_cys_pep = non_cys_pep.append(non_cys, ignore_index = True)
            # this removes any duplicate protein accessions to avoid reprocessing##
            quant_data = quant_data[(quant_data['Master Protein Accessions'] != x)]
    logger.debug('CysPep: {}'.format(cys_pep.shape))
    logger.debug('NonCysPep: {}'.format(non_cys_pep.shape))
    return two_unique_cys, cys_pep, non_cys_pep


def consensus(output_dict, col_heads):
    """
    Creates consensus dataframe from output dictionary, using the
    columns listed in col_heads to merge on

    Parameters:
    output_dict: dictionary
        containing treatmentID (key) mapped to output_dataframes (values)
    col_heads: list
        list of column headers to be used to merge all dataframes.
        Must exist in all dataframes.

    Returns:
    all_samples: DataFrame
        output dataframe containing all sample dataframes merged on
        the columns listed in col_heads

    """
    all_samples = pd.DataFrame()
    #generating empty columns in all_treats to use for merge
    for col in col_heads:
        all_samples[col] = ''

    #merging each dataframe from output dictionary according to columns
    for sampleID, dataframe in output_dict.items():
        #merge_ready = dataframe.set_index('Sequence', drop=True)
        all_samples = all_samples.merge(dataframe, how='outer', on=col_heads)

    return all_samples


def filter_NaNs(dataframe, filter_type='consecutive', threshold=1):
    """
    Filters rows from the consensus_df that contain too many missing values,
    according to total or consecutive mode

    Parameters:
        dataframe: pandas dataframe
            dataframe containing descriptive and numerical columns to be
            filtered. Numerical columns must contain floats.
        filter_type: string
            may be either consecutive or total. Consecutive filters out rows
            which have two consecutive NaNs, while total filters out rows
            that do not have at least threshold values
        threshold: int
            if filter_type is total, number of NaNs allowed in numeric columns
            if filter_type is consecutive, number of NaNs allowed consecutively

    Returns:
        dataframe_filtered: DataFrame
            dataframe containing rows that meet the filtering condition
    """
    #filter_type = filter_type.lower

    if filter_type == 'total':
        logger.info('Total filtering active')
        #to adjust threshold number according to what is required by dropna
        #dropna is the number of NaNs allowed
        threshold = dataframe.shape[1] - threshold
        #to take only peptides with less that the threshold number of NaNs
        dataframe_filtered = dataframe.dropna(thresh=threshold)

    elif filter_type == 'consecutive':
        logger.info('Consecutive filtering active')
        #to filter out more than one consecutive NaN for each peptides
        dataframe_filtered = dataframe.copy()
        threshold += 1

        for index, row in dataframe.iterrows():
            for x in range (0, (len(row)-1)):
                if type(row[x]) is float:
                    if pd.isna(row[x:x+threshold]).values.all():
                        dataframe_filtered.drop(index, axis=0, inplace=True)
                        break

    else:
        dataframe_filtered = pd.DataFrame()
        logger.debug("Filter type was not recognised. Please try 'consecutive' or 'total'.")

    return dataframe_filtered


def colour_column(df, col_test, colour_name):
    """Assigns colour to results according to -Log10(p-value) determined by t-test above threshold of 1.3

    Parameters
    ----------
    df : dataframe
        Description of parameter `df`.
    col_test : string
        Name of column to test for significance
    colour_name : string
        Name of column to store colour output

    Returns
    -------
    dataframe
        Returns entire input dataframe, with appended colour column
    """


    df[colour_name] = 'blue'
    df.loc[(df[col_test] > 1.3), [colour_name]] = 'red'

    return df
