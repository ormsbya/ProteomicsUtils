import numpy as np
import pandas as pd
from scipy import stats
import os
import logging
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
