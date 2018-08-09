
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


def main(input_path, output_path, sample_name, sample_type='whole_cell', replicate_threshold=0, simple=True, interactive=True, Bokeh_plot=True):

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
    # remove any proteins not seen in all replicates
    protein_AR_summary = DataWrangling.filter_NaNs(protein_AR_summary, filter_type='total', threshold=replicate_threshold)
    protein_AR_summary.reset_index(inplace=True, drop=True)
    logger.info(f"Protein AR: {protein_AR_summary.head(5)}")

    # Normalising each dataset to the Median Peptide Abundance
    protein_NormAR = protein_AR_summary.copy()
    for col in col_list:
            protein_NormAR[col] = protein_AR_summary[col]/median_list[col]
    logger.info(f"Protein abundances normalised to median peptide abundance: {protein_NormAR.head(5)}")


    # If IP sample, take the Log2 of each sample for t-tests
    if sample_type == "IP":
        logger.info(f'{sample_type} sample detected.')
        protein_Log2 = protein_NormAR.copy()
        for col in col_list:
                protein_Log2[col] = np.log2(protein_NormAR[col])
        logger.info(f"Log2 of Protein normalisaed abundances calculated: {protein_Log2.head(5)}")
    else:
        logger.info(f'{sample_type} sample detected. Using normalised abundances for one sample t-test')



    # Complete one-sample t-test on each row of NormProtAR using t-test_1samp function
    if sample_type == 'whole_cell':
        popmean = 1
        df = protein_NormAR
    elif sample_type == 'IP':
        popmean = 0
        df = protein_Log2

    logger.info(f"Calculating One Sample t-test with population mean {popmean}")
    df = CalcUtils.t_test_1samp(df, popmean, col_list)

    # Calculating the average abundance ratio
    logger.info(f"Calculating mean normalised Abundance Ratio...")
    df = CalcUtils.row_mean(df, col_list, 'Average')

    # Appending other columns of interest for the volcano plot
    # A volcano plot is constructed by plotting the negative log
    # of the p value on the y axis (usually base 10). This results
    # in data points with low p values (highly significant) appearing toward the top of the plot.
    df['Log10 p-val'] = -(np.log10(df['p-value']))
    if sample_type == 'whole_cell':
        logger.info(f"Calculating Log2 Average normalised Abundance Ratio, and -Log10(p-value)...")
        df['Log2 Av AR'] = np.log2(df['Average'])

    elif sample_type == 'IP':
        logger.info(f"Calculating Average Log2 normalised Abundance Ratio, and -Log10(p-value)...")
        df['Log2 Av AR'] = protein_Log2['Average']

    # To produce the colour column, change x and y limits in original function
    logger.info(f"Producing colour column...")
    xcol = 'Log2 Av AR'
    ycol = 'Log10 p-val'
    df = CalcUtils.colour_column_volc(df, xcol, ycol)
    logger.info(f"Post calculation results: {df.head(5)}")

    # Collecting dataframes and descriptors to save using the df_to_excel function
    data_frames = [calcs, protein_AR_summary, protein_NormAR, df]
    sheetnames = ['Med+Mean Calcs', 'Protein AR', 'ProtAR Norm to Med', 'Significance_test']
    output = output_path+sample_name+'ProteinAbundance_Results.xlsx'
    FileHandling.df_to_excel(output, sheetnames, data_frames)
    logger.info(f"Dataframes saved to excel file at {output}...")

    logger.info(f"Preparing data for volcano plot")
    # Gathering data for the scatter (volcano) plot
    xdata, xlabel = (df['Log2 Av AR'], 'Log2 Av. Abundance Ratio')
    ydata, ylabel = (df['Log10 p-val'], '-Log10 p-value')
    title = sample_name
    datalabels = df['Accession']
    colours = df['colours']


    if simple:
        # for simple volcano plot, which is saved into the pdf
        fig1 = PlotUtils.simple_scatter(xdata,ydata,title, xlabel, ylabel, colours)
        output = output_path+sample_name+'Simple_Volcano_'
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
        output = output_path+sample_name+'Interactive_Volcano_'
        FileHandling.fig_to_pdf([fig2], output)
        FileHandling.fig_to_svg(['Interactive_Volcano'],[fig2], output)


    if Bokeh_plot:
        output = output_path+sample_name+"_VolcanoPlot_Bokeh.html"
        output_file(output, title=sample_name)
        logger.info(f"Output html will be saved to {output_path}")

        hovers = [('Protein', '@Accession'),
            ('Gene', '@Description'),]

        fig3 = PlotUtils.bokeh_volcano_maker(df=df, c_col='Log10 p-val', y_col='Log10 p-val', x_col='Log2 Av AR', title=sample_name+' Volcano Plot', hover_list=hovers)
        show(fig3)

    # Saving figures to pdf and as svg files
    logger.info(f"Volcano plots saved to {output_path}")
    logger.info(f"Analysis complete for {sample_name}")


if __name__ == "__main__":
    #default parameters if no command line arguements given
    input_path = 'C:/Users/dezer_000/Documents/App_Dev_Projects/MS_Urea_Analysis_Project/TPEdenat/test_data/DC113-DC120_Compiled.xlsx'
    output_path = 'C:/Users/dezer_000/Documents/App_Dev_Projects/MS_Urea_Analysis_Project/TPEdenat/test_data/'
    sample_name = 'MG132'
    sample_type = 'whole_cell'
    replicate_threshold = 3
    main(input_path, output_path, sample_name, sample_type, replicate_threshold)
