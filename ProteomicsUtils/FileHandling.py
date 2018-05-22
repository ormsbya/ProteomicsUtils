"""
General collection of functions for handling input and output of data, including reading excel files, iterating over folders, saving dataframes to excel and saving figures to pdf or svg formats.
"""

import pandas as pd
import os
import re
from matplotlib.backends.backend_pdf import PdfPages
import xlrd
import logging
from ProteomicsUtils.LoggerConfig import logger_config

logger = logger_config(__name__)
logger.info('Import ok')

def dialog_path_finder():
    """ Antiquated function for opening a TKinter dialog window to receive user input selecting the input file. Splits file path to collect information on sample being processed.

    Returns
    -------
    input_path, output_path, sample_label: tuple
        tuple of strings containing input and output paths, and sample name.

    """
    #Create a reference to annoying parent window
    root = tkinter.Tk()
    #Close annoying parent window
    root.withdraw()
    #Get user input to select file
    input_path = filedialog.askopenfilename()
    results = input_path.split('/')
    sample_label = results[-1].split('.')
    sample_label = sample_label[0]
    logger.debug('Sample Name: {}'.format(sample_label))
    #logger.debug (results)
    resultpath = results[0:-1]
    output_path = '/'.join(resultpath) +'/'+sample_label+ '_'
#    logger.debug (output)

    return (input_path, output_path, sample_label)


def file_reader(path):
    """Antiquated function for opening excel files which are difficult (as are those that come from ProteomeDiscoverer!) and storing resulting data in DataFrame.

    Parameters
    ----------
    path : str
        full path to xlsx file to be opened.

    Returns
    -------
    dataframe
        containing data from first sheet provided by input path.

    """
    filename = path
    data_frame = pd.read_excel(
        filename,
        delimiter=u'\t', encoding='utf-8', skiprows=0,
        keep_default_na=False, na_values=['NA', 'N/A', 'nan', 'NaN', 'NULL', ''], comment=None,
        header=0, thousands=None, skipinitialspace=True,
        mangle_dupe_cols=True, quotechar='"',
        index_col=False)
    logger.debug('Total_data: {}'.format(data_frame.shape))
    return data_frame


def folder_iterator(input_folder, do_funcs, fileext = '.xlsx'):
    """
    Iterates through a folder, performs a set of tasks on each
    file according to file extension.

    Parameters:
    input_folder: string
        directory to iterate through, without final file seperator
    dofuncs: function
        contains functions to be applied to each file found by file iterator. Should be defined in __main__ and passed to folder iterator.
    fileext: string
        extension type of files to be processed in .ext format

    Returns:
    output_dict: dictionary
        mapped output from do_funcs to individual sample name
    """
    sample_list = []
    output_dict = {}

    for file in os.listdir(input_folder):
        if file.endswith(fileext):
            logger.info('Processing {}'.format(file))
            #grab sample_namename according to filename
            filename, extension = os.path.splitext(file)
            sample_name = filename.split('_')[-2]
            sample_list.append(sample_name)
            #create output results path using directory
            output_folder = input_folder+'/Results/'
            if not os.path.exists(output_folder):
                os.makedirs(output_folder)
            output_path = output_folder+sample_name+'_Results'+fileext
            input_path=(input_folder+'/'+file)
            #initiate empty list to store any problematic files
            exceptions = []
            try:
                #pass the file and output_path to the do_funcs function for processing
                output_per_treat = do_funcs(input_path, output_path, sample_name)
                #take return from do_funcs and append to dict
                output_dict[sample_name] = output_per_treat
                logger.info('Analysis completed for {}'.format(sample_name))
            except xlrd.XLRDError:
                logger.debug('File {} unable to be processed'.format(filename))
                pass

    return output_dict

def PD_compiler(input_folder, fileext='.xlsx'):
    """Creates compiled xlsx file from seperate Peptides and Proteins files exported from ProteomeDiscoverer.

    Parameters
    ----------
    input_folder : str
        Full path to folder containing xlsx files to be compiled.
    fileext : str (Default: 'xlsx')
        file extension of the files to be collected and compiled

    Returns
    -------
    output_dict
        dictionary containing filenames mapped to list of sample_name, result_type, dataframe for each file. Result type is determined from input file name.

    """
    sample_list = []
    output_dict = {}

    for file in os.listdir(input_folder):
        if file.endswith(fileext):
            logger.info('Processing {}'.format(file))
            input_path=(input_folder+'/'+file)
            #grab sample_name according to filename
            filename, extension = os.path.splitext(file)
            result_type = filename.split('_')[-1]
            sample_name = filename.split('_')[-2]
            logger.debug('Sample {}, Type {}'.format(sample_name, result_type))
            sample_list.append(sample_name)
            output_path = input_folder+sample_name+'_Compiled'+fileext

            #initiate empty list to store any problematic files
            exceptions = []
            try:
                #pass the file to the file_reader function for processing
                dataframe = file_reader(input_path)
                #take return from do_funcs and append to dict
                output_dict[filename] = [sample_name, result_type, dataframe]
                logger.info('{} dataframe collected'.format(filename))
            except xlrd.XLRDError:
                logger.debug('File {} unable to be processed'.format(filename))
                pass

    #iterate through unique sample names
    for sample in set(sample_list):
        #initialising empty lists
        sheetnames = []
        dataframes = []
        #collect all values from the dict that contain the sample name
        list_of_data = [value for key, value in output_dict.items() if sample in key]
        #for each value in the list, append the sample name, result type and data into lists
        for data in list_of_data:
            output_path = input_folder+'/'+data[0]+'_Compiled'+fileext
            sheetnames.append(data[1])
            dataframes.append(data[2])
        #Using utils function, save each sample to a single excel file with Pep and Prot labelled sheets
        df_to_excel(output_path, sheetnames, dataframes)
        logger.debug('Compiled {} saved to {}'.format(sample, output_path))
    return output_dict

def sheet_reader(path, sheetname):
    """Reads specific sheet from xlsx file located at path into DataFrame."""

    total_data = pd.read_excel(path, sheetname=sheetname)
    return total_data



def df_to_excel(output_path, sheetnames, data_frames):
    """Saves list of dataframes to a single excel (xlsx) file.

    Parameters
    ----------
    output_path : str
        Full path to which xlsx file will be saved.
    sheetnames : list of str
        descriptive list of dataframe content, used to label sheets in xlsx file.
    data_frames : list of DataFrames
        DataFrames to be saved. List order must match order of names provided in sheetname.

    Returns
    -------
    None.

    """
    if not output_path.endswith('.xlsx'):
        output_path = output_path+'Results.xlsx'
    writer = pd.ExcelWriter(output_path, engine='xlsxwriter')
    # Convert the dataframe to an XlsxWriter Excel object.
    for x in range(0, len(sheetnames)):
        sheetname = sheetnames[x]
        data_frame = data_frames[x]
        data_frame.to_excel(writer, sheet_name=sheetname, index=False)
    # Close the Pandas Excel writer and output the Excel file.
    writer.save()


def fig_to_pdf(figs, output_path, fig_type=None):
    """Save matplotlib figure objects to pdf document.

    Parameters
    ----------
    figs : dict or list
        Container with figure objects, either list or dictionary with figure objects as keys.
    output_path : str
        Partial path to folder to which figures will be saved. Figures.pdf is appended internally.
    fig_type : str, optional
        Appended to output_path prior to 'Figures.pdf' if provided.

    Returns
    -------
    None.

    """
    if fig_type:
        output_path = output_path+fig_type+'_'
    if isinstance(figs, dict):
        logger.info('Figure dictionary found')
        figs = list(figs.values())
    logger.info(figs)
    # page manager to allow saving multiple graphs to single pdf

    pdf = PdfPages(output_path + 'Figures.pdf')
    for fig in figs:
        pdf.savefig(fig)
        # fig.clf() #prevent plots from being overlayed onto the first
    # close pdfpages to allow open access
    pdf.close()


def fig_to_svg(fig_names, fig_list, output_path):
    """Save matplotlib figure objects to svg documents.

    Parameters
    ----------
    fig_names : list
        names given to the output svg files in the file path
    fig_list : list
        List of matplotlib figure objects, in order corresponding to fig_names.
    output_path : str
        Partial path to folder to which figures will be saved. Fig_name and extension (.svg) are appended internally.

    Returns
    -------
    None.

    """
    x = 0
    # figs = [plt.figure(n) for n in plt.get_fignums()]
    for fig in fig_list:
        figname = fig_names[x]
        filename = output_path + figname
        fig.savefig(filename + '.svg', transparent=True)
        x += 1
