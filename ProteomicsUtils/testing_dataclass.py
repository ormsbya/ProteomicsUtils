from dataclasses import dataclass
from typing import List
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from ProteomicsUtils.LoggerConfig import logger_config

logger = logger_config(__name__)
logger.info('Import successfull')

@dataclass
class Peptide:
    sequence: str
    abundance_ratio: List[float]
    accession: str

    def simple_scatter(self):
        x_vals = np.arange(0, len(self.abundance_ratio), 1)
        fig = plt.figure()
        plt.scatter(x_vals, self.abundance_ratio)
        return fig


@dataclass
class Protein:
    accession: str
    abundance_ratio: List[float]
    num_peptides: int
    peptides: List[Peptide]

test_data = pd.read_excel('C:/Users/dezer_000/Desktop/Current Analysis/180619_MQ_TPE+MG132/combined/Python_Analysis/MQ_clean_test_data_sample.xlsx')

test_data

peptide_1 = Peptide('QPNG', [0.1, 0.2], 'P14873')
peptide_2 = Peptide('ABCD', [0.3, 0.4], 'P15429')
peptide_3 = Peptide('EFGH', [0.5, 0.6], 'P14873')
peptide_4 = Peptide('IJKL', [0.7, 0.8], 'P15429')

protein_1 = Protein('P14873', [0.976, 0.251], 2, [peptide for peptide.accession == 'P14873'])

protein_2 = Protein('P15429', [0.967, 0.215], 2, [])


abundance_ratio_cols = [col for col in test_data.columns.tolist() if 'Replicate ' in col]
sequence = 'Sequence'
accession = 'ProteinID'

peptides = []
for i, peptide in test_data.iterrows():
    logger.info(f'Processing row {i} containing {peptide[sequence]}')
    peptides.append(Peptide(peptide[sequence], peptide[abundance_ratio_cols], peptide[accession]))
peptides


for peptide in peptides:
    logger.info(f'{peptide.accession}')

for peptide in peptides:
    fig = peptide.simple_scatter()
    plt.show(fig)
    
