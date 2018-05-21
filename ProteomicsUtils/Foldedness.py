import os, sys
import logging
import numpy as np
from LoggerConfig import logger_config
#Personal Modules
import StatUtils, CalcUtils, FileHandling, DataWrangling, PlotUtils
import AtomPandasConfig

logger = logger_config(__name__)
logger.info("Import Successful")

def volcano_plot(summary_data):
    xdata, xlabel = (summary_data['Log2 Average Ratio'], 'Log2 Average Cys/NonCys Ratio')
    ydata, ylabel = (summary_data['-Log10 p-Value'], '-Log10 p-Value')
    title = sample_name
    colours = summary_data['p-value colour']
    fig = PlotUtils.simple_scatter(xdata,ydata,title, xlabel, ylabel, colours)
    return fig

def foldedness_scatter(summary_data):
    xdata, xlabel = (summary_data.loc[(summary_data['p-value'] <0.05), ['Log2 Average NC']], 'Log2 NonCys Abundance')
    ydata, ylabel = (summary_data.loc[(summary_data['p-value'] <0.05), ['Log2 Average Ratio']], 'Log2 Cys/NonCys')
    title = sample_name
    colours = summary_data.loc[(summary_data['p-value'] <0.05), 'p-value colour']
    fig = PlotUtils.simple_scatter(xdata,ydata,title, xlabel, ylabel, colours)
    return fig

def do_funcs(input_path, output_path, sample_name, do_plots=True):
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
    logger.info('Input Folder: {}'.format(input_folder))
    logger.info('Input Path: {}'.format(input_path))

    total_data = FileHandling.file_reader(input_path)
    quant_data, col_list = DataWrangling.quantified_data(total_data)
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
    logger.info(ratio_col)
    NC_col = [col for col in abundance_cols if '_NC' in col]
    logger.info(NC_col)
    C_col =  [col for col in abundance_cols if ')_' not in col]
    logger.info(C_col)

    #collect only columns of interest for summary table
    select_col = ['Master Protein Accessions', 'Annotated Sequence'] + abundance_cols
    summary_data = summary_table[select_col]
    summary_data.dropna(axis=0, how='any', thresh=None, subset=abundance_cols, inplace=True)
    summary_data = summary_data.reset_index(drop=True)
    logger.info(summary_data)


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
    logger.info(summary_data)

    #####for peptides seen more than once in a sample, take average ratio to give only unique ratios for each peptide
    av_summary = CalcUtils.single_pep_av(summary_data)

    #Saving all dataframes so far to excel results document
    data_frames = [total_data, quant_data, two_unique_cys, cys_pep, non_cys_pep, summary_table, summary_data, av_summary]
    sheetnames = ['Total Data', 'Quant Data', 'TwoUniqueCYS', 'CysPep', 'NonCysPep', 'Summary Info', 'Summary Data', 'Single Peptide Average']
    FileHandling.df_to_excel(output_path, sheetnames, data_frames)
    logger.info("All dataframes saved to {}".format(output_path))

    if do_plots:
        pass


    return summary_data


def main(input_path):
    """Utility function to check the format of the input path, and applies functions as appropriate.
    Parameters
    ----------
    input_path : string
        Full path to either folder or file to be processed.
    Returns
    -------
    None. Saves the associated results via the do_funcs function.
    """
    output_path = input_path+'/Results/'
    if os.path.isdir(input_path):
        output_dict = FileHandling.folder_iterator(input_folder, do_funcs)
    else:
        filename, extension = os.path.splitext(file)
        sample_name = filename.split('_')[-2]
        do_funcs(input_path, output_path, sample_name)


if __name__ == "__main__":
    #check_input(input_path)
    pass


##########################From here############################################
## Notes on appropriate format for input and output strings would be useful to
## avoid ResultsResults file naming!
input_folder = 'C:/Users/dezer_000/Desktop/Current Analysis/180503_Stress_Compilation/Raw Results'
input_path = 'C:/Users/dezer_000/Desktop/Current Analysis/180503_Stress_Compilation/Raw Results/BarnaseEx4_170901_Dezerae_Sample2_Peptides.xlsx'
output_path = 'C:/Users/dezer_000/Desktop/Current Analysis/180503_Stress_Compilation/Foldedness Results/'
sample_name = 'BarnaseEx4'
##########################To here############################################


summary_data = do_funcs(input_path, output_path, sample_name)
main(input_path)


###http://newcoder.io/api/part-4/
###Checkout this to create the argparse function, and how to link into main function

figures = {}


figures['Foldedness Scatter'] = foldedness_scatter(summary_data)
figures['Volcano'] = volcano_plot(summary_data)
import matplotlib.pyplot as plt
for key, value in figures.items():
    print(key)
    print(value)
    plt.show(value)



#### Save figs to pdf
FileHandling.fig_to_pdf(figures, output_path, fig_type='Volcano+Foldedness')
logger.info(f"Figures saved to {output_path}")


"""
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--giantbomb-api-key', required=True,
                        help='API key provided by Giantbomb.com')
    parser.add_argument('--cpi-file',
                        default=os.path.join(os.path.dirname(__file__),
                                             'CPIAUCSL.txt'),
                        help='Path to file containing the CPI data (also acts'
                             ' as target file if the data has to be downloaded'
                             'first).')
    parser.add_argument('--cpi-data-url', default=CPI_DATA_URL,
                        help='URL which should be used as CPI data source')
    parser.add_argument('--debug', default=False, action='store_true',
                        help='Increases the output level.')
    parser.add_argument('--csv-file',
                        help='Path to CSV file which should contain the data'
                             'output')
    parser.add_argument('--plot-file',
                        help='Path to the PNG file which should contain the'
                             'data output')
    parser.add_argument('--limit', type=int,
                        help='Number of recent platforms to be considered')
    opts = parser.parse_args()
    if not (opts.plot_file or opts.csv_file):
        parser.error("You have to specify either a --csv-file or --plot-file!")
    return opts

    """
