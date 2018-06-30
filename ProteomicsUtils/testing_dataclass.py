from dataclasses import dataclass
from typing import List


@dataclass
class Peptide:
    sequence: str
    abundance_ratio: List[float]

@dataclass
class Protein:
    accession: str
    abundance_ratio: List[float]
    num_peptides: int
    peptides: List[Peptide]

peptide_1 = Peptide('QPNG', [0.1, 0.2])
peptide_2 = Peptide('ABCD', [0.3, 0.4])
peptide_3 = Peptide('EFGH', [0.5, 0.6])
peptide_4 = Peptide('IJKL', [0.7, 0.8])

protein_1 = Protein('P14873', [0.976, 0.251], 2, [peptide_1, peptide_2])

protein_2 = Protein('P15429', [0.967, 0.215], 2, [peptide_3, peptide_4])
