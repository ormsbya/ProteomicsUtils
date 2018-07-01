from dataclasses import dataclass, field
from typing import List
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from ProteomicsUtils.LoggerConfig import logger_config

logger = logger_config(__name__)
logger.info('Import successfull')


## Creating Peptide and Protein classes
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

    def contains_cys(self):
        cys = False
        if 'C' in self.sequence:
            cys = True
        return cys


@dataclass
class Protein:
    accession: str
    abundance_ratio: List[float]  = field(default_factory=list)
    peptides: List[Peptide]  = field(default_factory=list)

## Grabbing example cleaned test data
test_data = pd.read_excel('C:/Users/dezer_000/Desktop/Current Analysis/180619_MQ_TPE+MG132/combined/Python_Analysis/MQ_clean_test_data_sample.xlsx')

## Example instantiating Peptide class to specific variable
peptide_1 = Peptide('QPNG', [0.1, 0.2], 'P14873')


## Collecting column descriptors for different properties
abundance_ratio_cols = [col for col in test_data.columns.tolist() if 'Replicate ' in col]
sequence = 'Sequence'
accession = 'ProteinID'


## Instantiating peptides in a list as instances of the Peptides class
peptides = []
for i, peptide in test_data.iterrows():
    logger.info(f'Processing row {i} containing {peptide[sequence]}')
    peptides.append(Peptide(peptide[sequence], peptide[abundance_ratio_cols], peptide[accession]))

## To grab the accession number for each peptide
for peptide in peptides:
    logger.info(f'{peptide.accession}')

## E.g. using a class method to create individual plots for each peptide
for peptide in peptides:
    fig = peptide.simple_scatter()
    plt.show(fig)

## Example collecting info for a single protein
protein_1 = Protein('Q80X82')
protein_1.peptides = [peptide.sequence for peptide in peptides if peptide.accession == 'Q80X82']

## To instantiate proteins in a list
proteins = []
for protein in test_data[accession].unique():
    logger.info(f'Processing {protein}')
    proteins.append(Protein(protein, [], [peptide for peptide in peptides if peptide.accession == protein]))

## Testing ease of collecting specific information for one protein
protein_of_interest = [protein for protein in proteins if protein.accession == 'Q80X82'][0]
for peptide in protein_of_interest.peptides:
    print (peptide.sequence)

## checking if peptides contain a cys residue

for peptide in peptides:
    cys = peptide.contains_cys()
    if cys:
        logger.info(f'Cys found in {peptide.sequence}')
    else:
        logger.info(f'No cys found in {peptide.sequence}')
