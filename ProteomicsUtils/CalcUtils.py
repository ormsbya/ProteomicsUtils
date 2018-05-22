"""
General collection of functions for calculating .
"""

import numpy as np
import pandas as pd
from scipy import stats
import math
import os
import logging
from ProteomicsUtils.LoggerConfig import logger_config

logger = logger_config(__name__)
logger.info('Import ok')


def cys_div_noncys(cys_pep, non_cys_Av, col_list):
    """
    Calculates the cys/noncys ratio for each cysteine peptide using
    non_cys_Av per protein values.

    Parameters:
    cys_pep: dataframe containing at minimum ProteinID and Abundance Ratio
        columns for each cysteine-containing peptide
    non_cys_Av: dataframe containing ProteinID and average non-cysteine
        peptide abundance ratios
    col_list:   list containing headers for all columns containing
        abundance ratios to be operated on

    Returns:
    cys_pep: df containing two new columns with non_cys_Av and
        Cys/NonCys ratio for each peptide

    """
    for x in col_list:
        #create replicate-specific col headers
        non_cys_header = str(x)+'_NC'
        cys_div_nc_header = str(x)+'_Cys/NonCys'
        #append non-cys ratio
        for index, Series in cys_pep.iterrows():
            cys_pep.loc[index, non_cys_header] = non_cys_Av.loc[index, x]
        #to calculate cys/non-cys
        cys_pep[cys_div_nc_header] = cys_pep[x]/cys_pep[non_cys_header]

    return cys_pep

def single_element_av(dataframe, col_name):
    """
    Creates summary dataframe where each peptide has a single
    ratio value (if seen multiple times in a single sample),
    by averaging the ratio for peptides observed multiple times

    Parameters:
    dataframe: DataFrame
        summary dataframe containing ProteinID, Sequence and Ratio
        columns

    Returns:
    av_summary: DataFrame
        summary dataframe where each peptide appears only one.
    """
    #collect unique peptides
    unique_index = dataframe[col_name].unique()
    logger.info(f"Total unique elements found: {len(unique_index)}")
    #create empty dataframe
    av_summary = pd.DataFrame()
    #get column names containing Ratio
    col_heads = [col for col in dataframe.columns if '_Cys/NonCys' in col]
    logger.info(f"Columns for calculation: {col_heads}")

    for element in unique_index:
        #collect all instances of a single peptide
        ratios = dataframe.loc[(dataframe[col_name] == element)]
        #check for 1 entry, or multiple entries
        if len(ratios) == 1:
            #if only one ratio was found for a peptide, take that
            av_summary = av_summary.append(ratios, ignore_index = True)
        if len(ratios) > 1:
            #if more than one ratio was found for a peptide, take average
            for column in col_heads:
                av = ratios.iloc[0]
                av[column] = ratios[column].mean()
            av_summary = av_summary.append(av)

    return av_summary


def non_cys_AR(cys_summ, non_cys_summ):
    """Calculates the average abundance ratio of non-cysteine containing peptides from a given protein.

    Parameters
    ----------
    cys_summ, non_cys_summ : DataFrames
        DataFrames containing data for the (1) cysteine, and (2) non-cysteine peptides respectively.

    Returns
    -------
    DataFrame
        New dataframe with protein means appended in abundance ratio columns for each replicate.

    """

    non_cys_means = pd.DataFrame()
    for x in non_cys_summ['Master Protein Accessions']:
        protein = non_cys_summ.loc[(non_cys_summ['Master Protein Accessions'] == x)]
        col_heads = [col for col in protein.columns if 'Abundance Ratio: (' in col]

        for l in range (0, len(col_heads)):
            column = col_heads[l]
            if protein.shape[0] == 1:
                non_cys_means.loc[x,column] = non_cys_summ.loc[x, column]
            if protein.shape[0] > 1:
                #logger.debug(protein[column])
                non_cys_mean = protein[column].mean()
                non_cys_means.loc[x,column] = non_cys_mean
            #NonCysAbun_dict[x] = rep_means
        #this removes any duplicate protein accessions to avoid reprocessing##
        non_cys_summ = non_cys_summ[(non_cys_summ['Master Protein Accessions'] != x)]
    #logging.debug (non_cys_means)
    return non_cys_means


def colour_column_volc(df, xcol, ycol):
    """Generates red/blue/gray colour column for volcano plot data according to the value in x and y cols for each protein."""
    df['colours'] = 'gray'
    df.loc[(df[xcol] > 1) & (df[ycol] > 1.3), ['colours']] = 'red'
    df.loc[(df[xcol] < -1) & (df[ycol] > 1.3), ['colours']] = 'blue'
    return df



def mean_med_calc(vals, samplename):
    """Calculates the mean and median for a seet of values, appended to dictionary mapped to samplename"""
    vals = vals.dropna()
    mean_val = np.mean(vals)
    median_val = np.median(vals)
    calcs_dict[samplename] = [mean_val, median_val]
    return calcs_dict


# def row_mean_name(df, colstart, colstop, colname):
#     colstart = df.columns.get_loc(colstart)
#     colstop = df.columns.get_loc(colstop)
#     for i in range(df.shape[0]):
#         vals = df.iloc[i, colstart:(colstop + 1)]
#         vals = vals.dropna()
#         mean_val = np.mean(vals)
#         df.loc[i, colname] = mean_val
#     return df


def cys_abun_change(cys_pep, non_cys_pep):
    """Calculates the average abundance ratio (and log2) of cysteine containing peptides normalised to non-cysteine containing peptides from that protein.

    Parameters
    ----------
    cys_pep, non_cys_pep : DataFrames
        DataFrames containing data for the (1) cysteine, and (2) non-cysteine peptides respectively.

    Returns
    -------
    (AbundanceRatios, Log2Ratios, cys_pep) : Tuple
        New dataframe with cys peptide AR and Cys/Non-cys AR, plus Log2 of ratios, appended in new DataFrames.
    """

    CysAbun_dict = {}
    NonCysAbun_dict = {}
    Change_dict = {}
    changelist = []
    non_cys_list = []
    cys_pep_per_protein = cys_pep
    cys_pep_per_peptide = cys_pep

    for x in cys_pep_per_protein['Master Protein Accessions']:
        protein = cys_pep_per_protein.loc[(cys_pep_per_protein['Master Protein Accessions'] == x)]
        if protein.empty != True:
            #to process a multi-consensus file (in which the av. abundance has been calc)
            col_head = [col for col in protein.columns if 'Abundance Ratio (Average)' in col]
            col_head = str(col_head[0])
            cys_mean = protein[col_head].mean()
            CysAbun_dict[x] = cys_mean
        #this removes any duplicate protein accessions to avoid reprocessing##
        cys_pep_per_protein = cys_pep_per_protein[(cys_pep_per_protein['Master Protein Accessions'] != x)]

    for x in non_cys_pep['Master Protein Accessions']:
        protein = non_cys_pep.loc[(non_cys_pep['Master Protein Accessions'] == x)]
        if protein.empty != True:
            non_cys_mean = protein[col_head].mean()
            NonCysAbun_dict[x] = non_cys_mean
        #this removes any duplicate protein accessions to avoid reprocessing##
        non_cys_pep = non_cys_pep[(non_cys_pep['Master Protein Accessions'] != x)]
    print ('Column used for calculations: ', col_head)

    for index, row in cys_pep_per_peptide.iterrows():
        proteinID = row['Master Protein Accessions']
        non_cys_av = NonCysAbun_dict[proteinID]
        change_ratio = row[col_head]/non_cys_av
        non_cys_list.append(non_cys_av)
        changelist.append(change_ratio)

    for x in CysAbun_dict:
        Protchange = CysAbun_dict[x] / NonCysAbun_dict[x]
        Change_dict[x] = Protchange

    ChangeRatios =  pd.DataFrame.from_dict(Change_dict, orient='index')
    ChangeRatios.columns = ['Cys/NonCys']
    CysRatios =  pd.DataFrame.from_dict(CysAbun_dict, orient='index')
    CysRatios.columns = ['Cys']
    NonCysRatios =  pd.DataFrame.from_dict(NonCysAbun_dict, orient='index')
    NonCysRatios.columns = ['NonCys']

    #add cys change per peptide and log2 of cys change per pep to cys_pep dataframe
    cys_pep['NonCys Av'] = non_cys_list
    cys_pep['Log2 NonCys Av'] = np.log2(non_cys_list)
    cys_pep['Change per peptide'] = changelist
    cys_pep['Log2 Change per peptide'] = np.log2(changelist)

    #joining abundance ratios per protein
    AbundanceRatios = CysRatios.join(NonCysRatios, how = 'outer')
    AbundanceRatios = AbundanceRatios.join(ChangeRatios, how = 'outer')

    #calculate mean abundance ratios
    cys_mean = AbundanceRatios['Cys'].mean()
    noncys_mean = AbundanceRatios['NonCys'].mean()
    change_mean = AbundanceRatios['Cys/NonCys'].mean()
    AbundanceRatios.loc['Mean'] = (cys_mean, noncys_mean, change_mean)

    #Calculate median abundance ratios
    cys_median = AbundanceRatios['Cys'].median()
    noncys_median = AbundanceRatios['NonCys'].median()
    change_median = AbundanceRatios['Cys/NonCys'].median()
    AbundanceRatios.loc['Median'] = (cys_median, noncys_median,change_median)

     #normalise NonCys to the overall median of NonCys peptides
    norm_median = AbundanceRatios['NonCys']/noncys_median
    norm_median = pd.DataFrame(norm_median)
    norm_median.columns = ['NonCys Median Norm']
    AbundanceRatios = AbundanceRatios.join(norm_median, how = 'outer')

    #to take log2 of all calculated abundance ratios
    Log2Ratios = np.log2(AbundanceRatios)

    return (AbundanceRatios, Log2Ratios, cys_pep)


def t_test_pair(df, cols_a, cols_b):
    """Completes paired t-test for each row of values, comparing those in col_a to those in col_b
    Parameters
    ----------
    df : dataframe
        Contains the complete dataset of interest
    cols_a : list
        list of column names to be included in the first set of values
    cols_b : list
        list of column names to be included in the second set of values
    Returns
    -------
    dataframe
        Returns entire input dataframe, with appended t-stat and p-value columns
    """

    for i in range(df.shape[0]):
        vals_a = df.loc[i, cols_a]
        vals_b = df.loc[i, cols_b]

        t_test_vals = stats.ttest_rel(vals_a, vals_b, axis=0)
        df.loc[i, 't-stat'] = t_test_vals[0]
        df.loc[i, 'p-value'] = t_test_vals[1]

    return df


def t_test_1samp(df, popmean, cols):
    """Completes one-sample t-test for each row of values, comparing those in cols to the popmean
    Parameters
    ----------
    df : dataframe
        Contains the complete dataset of interest
    cols : list
        list of column names to be included in the values
    popmean : float
        value to compare to (i.e. should be 0 or 1 for most ratio queries according to whether data has been logged (0) or not (1) prior to test)
    Returns
    -------
    dataframe
        Returns entire input dataframe, with appended t-stat and p-value columns
    """
    for i in range(df.shape[0]):
        vals = df.loc[i, cols]
        vals = vals.dropna()
        t_test_vals = stats.ttest_1samp(vals, popmean)
        df.loc[i, 't-stat'] = t_test_vals[0]
        df.loc[i, 'p-value'] = t_test_vals[1]
    return df

def row_mean(df, cols, result_name):
    """calculates mean of columns per row of a dataframe.
    Parameters
    ----------
    df : dataframe
        Description of parameter `df`.
    cols : list
        List of column names to be used for calculation
    result_name : list
        Name of column in which output is stored
    Returns
    -------
    dataframe
        Returns entire input dataframe, with appended mean columns
    """

    for i in range(df.shape[0]):
        vals = df.loc[i, cols]
        vals = vals.dropna()
        mean_val = np.mean(vals)
        df.loc[i, result_name] = mean_val
    return df
