
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit
#Personal Modules
from ProteomicsUtils import FileHandling, StatUtils
from ProteomicsUtils.LoggerConfig import logger_config


logger = logger_config(__name__)
logger.info('Import ok')

def per_protein_fitter(df_for_fitting, x_vals, group_col=None, element_col=None, svg=None):
    calculation_store = {}
    nonconverged = []
    if group_col:
        df_for_fitting = df_for_fitting.drop(group_col, axis=1)
        df_for_fitting.set_index(element_col, inplace=True)
    #to produce default x-values for fitting if none given
    if not x_vals:
        x_vals = np.arange(0, df_for_fitting.shape[1], 1)
        logger.info(f'No input for x_vals detected. Using {x_vals} instead.')
    for element in df_for_fitting.index.tolist():
        y_data = df_for_fitting.loc[element].tolist()
        logger.debug(f'x_vals: {x_vals}')
        logger.debug(f'y_vals: {y_data}')
        try:
            calculation_store[element] = StatUtils.sigmoid_calculator(x_vals, y_data)
        except(RuntimeError):
            nonconverged.append(element)
    fig = StatUtils.sigmoid_plotter(calculation_store)
    return fig

def main(input, output_path, sample_name, x_vals=None, test_element=None, group_col=None, element_col=None):
    pass

    logger.info(f"Analysing: {sample_name}")
    if not os.path.isdir(output_path):
        os.mkdir(output_path)
        logger.info(f"Output directory made at {output_path}")
    if isinstance(input, pd.DataFrame):
        logger.info('DataFrame input detected')
        raw_data = input
    elif os.path.isfile(input):
        logger.info(f'Input file being loaded from {input}')
        raw_data = pd.read_excel(input)
    else:
        logger.info(f"Incorrect input format detected. Please pass full file path, or a dataframe as input.")

    # Collects index column for plot names by default
    if not element_col:
        element_col = raw_data.index.tolist()
        logger.info(f'No element column detected. Using index instead.')


    #If a grouping column is passed (e.g. proteins) then all elements for that group are plotted together. Otherwise, individual elements are plotted separately, and are expected the be the index
    if group_col:
        col_for_dict = group_col
    else:
        col_for_dict = element_col

    per_element_dict = {}
    for element in raw_data[col_for_dict].unique():
        per_element_dict[element] = raw_data[raw_data[col_for_dict] == element]
    #If a test element is passed to function, then only produce the plot for that element, otherwise fit whole df
    if test_elements:
        data_for_fitting = {k: per_element_dict[k] for k in test_elements}
    else:
        data_for_fitting = per_element_dict
    #Create the figure for each group or element
    figure_dict = {}
    for key, value in per_element_dict.items():
        fig = per_protein_fitter(value, x_vals, group_col, element_col)
        fig.suptitle(key)
        figure_dict[key] = fig
        plt.show(fig)

    FileHandling.fig_to_pdf(figure_dict, output_path=output_path+sample_name, fig_type='_Sigmoids')
    if svg:
        FileHandling.fig_to_svg(fig_names=list(figure_dict.keys()), fig_list=list(figure_dict.values()), output_path=output_path+sample_name)
    logger.info(f'Figures saved to {output_path}')

    return figure_dict

if __name__ == "__main__":
    input_path = "C:/Users/dezer_000/Desktop/Current Analysis/180501_Urea_Exp8_Analysis/Thresholded_Results.xlsx"
    output_path = 'C:/Users/dezer_000/Desktop/'
    sample_name = 'Test_data'
    group_col = 'ProteinID'
    element_col = 'Sequence'
    main(input_path, output_path, sample_name, group_col=group_col, element_col=element_col)
