import pandas as pd
import seaborn as sns
import logging
from ProteomicsUtils.LoggerConfig import logger_config
import warnings

logger = logger_config(__name__)
logger.info('Import ok')


# setting the colour scheme to seaborn for plots
sns.set()
# forcing pandas to display df without wrapping lines and prevent raising annoying warning
pd.set_option('expand_frame_repr', False)
pd.options.mode.chained_assignment = None  # default='warn'

#Ignore FutureWarnings raised due to sheetnames
warnings.simplefilter(action='ignore', category=FutureWarning)

logger.info('Pandas display config complete')
