import os, sys
import logging
import numpy as np
import matplotlib.pyplot as plt
#Personal Modules
from ProteomicsUtils import StatUtils, CalcUtils, FileHandling, DataWrangling, PlotUtils
from ProteomicsUtils.LoggerConfig import logger_config

logger = logger_config(__name__)
logger.info("Import Successful")


def foldedness_scatter(summary_data):
    xdata, xlabel = (summary_data.loc[(summary_data['p-value'] <0.05), ['Log2 Average NC']], 'Log2 NonCys Abundance')
    ydata, ylabel = (summary_data.loc[(summary_data['p-value'] <0.05), ['Log2 Average Ratio']], 'Log2 Cys/NonCys')
    title = sample_name
    colours = summary_data.loc[(summary_data['p-value'] <0.05), 'p-value colour']
    fig = PlotUtils.simple_scatter(xdata,ydata,title, xlabel, ylabel, colours)
    return fig

def main(input_path, output_path, sample_name, do_plots=True):
    """
    Master function to apply a list of functions to the input file

    Parameters:
    input_path: string
        input path for the file to be processed
    output_path: string
        output path for which any output generated by functions will be saved
    sample_name: string
        sample name associated with the file to be processed.

    Returns:
    summary_table: DataFrame
        dataframe containing the summarised output of the functions
        applied in order
    """
    logger.info(f'Preparing to process: {sample_name}')
    logger.info(f"Input Path: {input_path}")

    total_data = FileHandling.file_reader(input_path)

    quant_data, col_list = DataWrangling.quantified_data(total_data)
    #raj only considers peptides that are observed in all replicates
    quant_data = quant_data.dropna(axis=0, how='any', thresh=None, subset=col_list)
    two_unique_cys, cys_pep, non_cys_pep = DataWrangling.Unique_Cys_sorter(quant_data)
    #set index of summary dataframes to the protein accession
    cys_pep = cys_pep.set_index(["Master Protein Accessions"], drop=False)
    non_cys_pep = non_cys_pep.set_index(["Master Protein Accessions"], drop=False)

    non_cys_Av = CalcUtils.non_cys_AR(cys_pep, non_cys_pep)

    summary_table = CalcUtils.cys_div_noncys(cys_pep, non_cys_Av, col_list)

    #collect list of columns for each type of ratio
    summary_table.reset_index(drop=True, inplace=True)
    abundance_cols =  [col for col in summary_table.columns if 'Abundance Ratio: (' in col]
    ratio_col = [col for col in abundance_cols if '_Cys/NonCys' in col]
    logger.info(f"Collecting ratio columns: {ratio_col}")
    NC_col = [col for col in abundance_cols if '_NC' in col]
    logger.info(f"Collecting NonCys columns: {NC_col}")
    C_col =  [col for col in abundance_cols if ')_' not in col]
    logger.info(f"Collecting Cys columns: {C_col}")

    #collect only columns of interest for summary table
    select_col = ['Master Protein Accessions', 'Annotated Sequence'] + abundance_cols
    summary_data = summary_table[select_col]
    summary_data.dropna(axis=0, how='any', thresh=None, subset=abundance_cols, inplace=True)
    summary_data = summary_data.reset_index(drop=True)
    logger.debug(F"Summary data acquired: {summary_data}")


    #rename columns to simple names
    summary_data = summary_data.rename(columns = {'Master Protein Accessions':'ProteinID',   'Annotated Sequence':'Sequence'})

    #####Paired T-test
    logger.info("Calculating t-test statistics")
    summary_data = CalcUtils.t_test_pair(summary_data, C_col, NC_col)

    #####Include -log10 of p-value
    summary_data['-Log10 p-Value'] = -np.log10(summary_data['p-value'])

    #####Include Log2 of Average ratios
    summary_data = CalcUtils.row_mean(summary_data, ratio_col, 'Av. Ratio')
    summary_data['Log2 Average Ratio'] = np.log2(summary_data['Av. Ratio'])

    ##### Include average col for NC
    summary_data = CalcUtils.row_mean(summary_data, NC_col, 'NC Average')
    summary_data['Log2 Average NC'] = np.log2(summary_data['NC Average'])

    #####Assign colour column
    summary_data = DataWrangling.colour_column(summary_data, '-Log10 p-Value', 'p-value colour')
    logger.info("Completed calculations with summary data table")
    logger.debug(f"Summary data table post calculations: {summary_data}")

    #####for peptides seen more than once in a sample, take average ratio to give only unique ratios for each peptide
    av_summary = CalcUtils.single_pep_av(summary_data)

    #Saving all dataframes so far to excel results document
    data_frames = [total_data, quant_data, two_unique_cys, cys_pep, non_cys_pep, summary_table, summary_data, av_summary]
    sheetnames = ['Total Data', 'Quant Data', 'TwoUniqueCYS', 'CysPep', 'NonCysPep', 'Summary Info', 'Summary Data', 'Single Peptide Average']
    FileHandling.df_to_excel(output_path=output_path+sample_name, sheetnames=sheetnames, data_frames=data_frames)
    logger.info("All dataframes saved to {output_path}")

    if do_plots:
        logger.info(f"Preparing foldedness scatterplot for {sample_name}")
        figures = {}
        figures['Foldedness Scatter'] = foldedness_scatter(summary_data)
        #### Save figs to pdf
        FileHandling.fig_to_pdf(figures, output_path=output_path, fig_type=sample_name+'_Foldedness')
        logger.info(f"Figures saved to {output_path}")
        for key, value in figures.items():
            plt.show(value)

    return summary_data



if __name__ == "__main__":
    #default path to test data
    input_path = 'C:/Users/dezer_000/Documents/App_Dev_Projects/ProteomicsUtils/test_data/test_data_Peptides.xlsx'
    output_path = 'C:/Users/dezer_000/Documents/App_Dev_Projects/ProteomicsUtils/test_data/'
    sample_name = 'test_data'
    main(input_path, output_path, sample_name, do_plots=True)
