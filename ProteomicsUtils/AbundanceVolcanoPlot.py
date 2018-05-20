
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from bokeh.plotting import figure, show, output_file
from ProteomicsUtils import FileHandling, DataWrangling, PlotUtils, CalcUtils
from ProteomicsUtils.LoggerConfig import logger_config

logger = logger_config(__name__)
logger.info('Import Successful')


# setting the colour scheme to seaborn for plots
sns.set()


def main(input_path, output_path, sample_name, simple=True, interactive=True, Bokeh_plot=True):

    logger.info(f"Analysing: {sample_name}")
    if not os.path.isdir(output_path):
        os.mkdir(output_path)

    ## COLLECTING PEPTIDE ABUNDANCES FOR NORMALISATION ##
    # opening sheet by calling sheet_reader on path
    sheetname = 'Peptides'
    peptides_raw = pd.read_excel(input_path, sheetname)
    logger.info(f"Peptide data collected from {input_path}")
    logger.info(f"Number of peptides detected: {peptides_raw.shape[0]}")

    #collect list of columns containing abundance Ratios
    col_list = [col for col in peptides_raw.columns if 'Abundance Ratio: (' in col]
    logger.info(f"Columns detected for analysis: {col_list}")

    # calculating mean and median for pep abundance using mean_median_calc function
    calcs_dict = {}
    for column in col_list:
        vals = peptides_raw[column].dropna()
        mean_val = np.mean(vals)
        median_val = np.median(vals)
        calcs_dict[column] = [mean_val, median_val]
        logger.info(f"Med. and Mean calculated for {column}")

    # convert calcs to dataframe and set column labels
    calcs = pd.DataFrame.from_dict(calcs_dict, orient='index')
    calcs.columns = ['Mean', 'Median']
    calcs = calcs.sort_index()
    median_list = calcs['Median']
    mean_list = calcs['Mean']

    logger.info(f"Median: {median_list}")
    logger.info(f"Mean: {mean_list}")

    ## COLLECTING PROTEIN ABUNDANCES FOR VOLCANO PLOT ##
    # opening all Protein excel sheets by calling sheet_reader on path
    sheetname = 'Proteins'
    proteins_raw = pd.read_excel(input_path, sheetname)
    logger.info(f"Proteins imported from {input_path}")

    # creating summary data_frame of original protein data
    summary_cols = ['Accession','Description'] + col_list
    logger.info(f"Columns for summary: {summary_cols}")
    protein_AR_summary = proteins_raw[summary_cols]
    logger.info(f"Protein AR: {protein_AR_summary.head(5)}")

    # Normalising each dataset to the Median Peptide Abundance
    protein_NormAR = protein_AR_summary.copy()
    for col in col_list:
            protein_NormAR[col] = protein_AR_summary[col]/median_list[col]
    logger.info(f"Protein abundances normalised to median peptide abundance: {protein_NormAR.head(5)}")


    # Complete one-sample t-test on each row of NormProtAR using t-test_1samp function
    popmean = 1
    logger.info(f"Calculating One Sample t-test with population mean {popmean}")
    protein_NormAR = CalcUtils.t_test_1samp(protein_NormAR, popmean, col_list)

    # Calculating the average abundance ratio
    logger.info(f"Calculating mean normalised Abundance Ratio...")
    protein_NormAR = CalcUtils.row_mean(protein_NormAR, col_list, 'Average')


    # Appending other columns of interest for the volcano plot
    # A volcano plot is constructed by plotting the negative log
    # of the p value on the y axis (usually base 10). This results
    # in data points with low p values (highly significant) appearing toward the top of the plot.
    logger.info(f"Calculating Log2 Average normalised Abundance Ratio, and -Log10(p-value)...")
    protein_NormAR['Log10 p-val'] = -(np.log10(protein_NormAR['p-value']))
    protein_NormAR['Log2 Av AR'] = np.log2(protein_NormAR['Average'])

    # To produce the colour column, change x and y limits in original function
    logger.info(f"Producing colour column...")
    xcol = 'Log2 Av AR'
    ycol = 'Log10 p-val'
    protein_NormAR = CalcUtils.colour_column_volc(protein_NormAR, xcol, ycol)
    logger.info(f"Post calculation results: {protein_NormAR.head(5)}")

    # Collecting dataframes and descriptors to save using the df_to_excel function
    data_frames = [calcs, protein_AR_summary, protein_NormAR]
    sheetnames = ['Med+Mean Calcs', 'Protein AR', 'ProtAR Norm to Med']
    output = output_path+'ProteinAbundance_Results.xlsx'
    FileHandling.df_to_excel(output, sheetnames, data_frames)
    logger.info(f"Dataframes saved to excel file at {output}...")

    logger.info(f"Preparing data for volcano plot")
    # Gathering data for the scatter (volcano) plot
    xdata, xlabel = (protein_NormAR['Log2 Av AR'], 'Log2 Av. Abundance Ratio')
    ydata, ylabel = (protein_NormAR['Log10 p-val'], '-Log10 p-value')
    title = sample_name
    datalabels = protein_NormAR['Accession']
    colours = protein_NormAR['colours']


    if simple:
        # for simple volcano plot, which is saved into the pdf
        fig1 = PlotUtils.simple_scatter(xdata,ydata,title, xlabel, ylabel, colours)
        output = output_path+'Simple_Volcano_'
        FileHandling.fig_to_pdf([fig1], output)
        FileHandling.fig_to_svg(['Simple_Volcano'],[fig1], output)
        plt.show(fig1)

    if interactive:
        # for interactive volcano plot
        # create the scatterplot
        fig2 = PlotUtils.inter_scatter(xdata, ydata, xlabel, ylabel, colours, title, datalabels)
        # initial drawing of the scatterplot
        plt.plot()
        logger.info("Interactive scatterplot done")
        # present the scatterplot
        #plt.show()
        output = output_path+'Interactive_Volcano_'
        FileHandling.fig_to_pdf([fig2], output)
        FileHandling.fig_to_svg(['Interactive_Volcano'],[fig2], output)


    if Bokeh_plot:
        output = output_path+"_VolcanoPlot_Bokeh.html"
        output_file(output, title=sample_name)
        logger.info(f"Output html will be saved to {output_path}")

        hovers = [('Protein', '@Accession'),
            ('Gene', '@Description'),]

        fig3 = PlotUtils.bokeh_volcano_maker(df=protein_NormAR, c_col='Log10 p-val', y_col='Log10 p-val', x_col='Log2 Av AR', title=sample_name+' Volcano Plot', hover_list=hovers)
        show(fig3)

    # Saving figures to pdf and as svg files
    logger.info(f"Volcano plots saved to {output_path}")
    logger.info(f"Analysis complete for {sample_name}")


if __name__ == "__main__":
    #default parameters if no command line arguements given
    input_path = 'C:/Users/dezer_000/Documents/App_Dev_Projects/MS_Urea_Analysis_Project/TPEdenat/test_data/DC113-DC120_Compiled.xlsx'
    output_path = 'C:/Users/dezer_000/Documents/App_Dev_Projects/MS_Urea_Analysis_Project/TPEdenat/test_data/'
    sample_name = 'MG132'
    main(input_path, output_path, sample_name)
