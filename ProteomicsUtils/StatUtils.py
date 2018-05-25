import numpy as np
import pandas as pd
from scipy import stats
import os
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from ProteomicsUtils.LoggerConfig import logger_config


logger = logger_config(__name__)
logger.info('Import ok')




def t_test_1samp(df, popmean, colstart, colstop):
    for i in range(df.shape[0]):
        vals = df.iloc[i, colstart:(colstop + 1)]
        vals = vals.dropna()
        t_test_vals = stats.ttest_1samp(vals, popmean)
        df.loc[i, 't-stat'] = t_test_vals[0]
        df.loc[i, 'p-value'] = t_test_vals[1]
    return df

from scipy.optimize import curve_fit

def sigmoid(x, x0, k, a, c):
     y = a / (1 + np.exp(-k*(x-x0))) + c
     return y

def sigmoid_calculator(xdata, ydata):
    popt, pcov = curve_fit(sigmoid, xdata, ydata)
    #print (popt)

    x = np.linspace(min(xdata), max(xdata), 100)
    y = sigmoid(x, *popt)
    return (x, y, xdata, ydata)

def sigmoid_plotter(sigmoid_dictionary):
    fig = plt.figure()
    for key, value in sigmoid_dictionary.items():
        x, y, x_data, y_data = value
        plt.scatter(x_data, y_data, label=key)
        plt.plot(x,y, label=key+'_fit')
        plt.xlim(min(x_data), max(x_data))
        plt.xlabel('Urea Concentration')
        plt.ylabel('Cys/NonCys Abundance Ratio')
    plt.legend(loc='best')
    return fig

def per_protein_fitter(protein_df):
    calculation_store = {}
    nonconverge_pep = []
    protein_df = protein_df.drop('ProteinID', axis=1).set_index('Sequence')
    for peptide in protein_df.index.tolist():
        y_data = protein_df.loc[peptide].tolist()
        try:
            calculation_store[peptide] = sigmoid_calculator(urea_conc, y_data)
        except(RuntimeError):
            nonconverge_pep.append(peptide)
    fig = sigmoid_plotter(calculation_store)
    return fig
#calculation_store
