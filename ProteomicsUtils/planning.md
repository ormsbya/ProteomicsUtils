## Thinking about using dataclass
 - Proteins would be the top level class
 - peptides would be a subclass that inherit all virtues of the parent class?

#### proteins
Available from PD data:
- Accession (str)
- abundance_ratio (list? of floats)
- num_peptides
- peptides (list of peptides from peptide data - optional)

Values that can be added by calculations:
- av_abundance_ratio
- log10_pvalue
- gene_name
- gene_ontology

Methods:
- average non-cys

#### peptides
- abundance_ratio (list? of floats)

Methods:
- contains cys
-


## Example implementation from [Real Python](https://realpython.com/python-data-classes/)

from dataclasses import dataclass
from typing import List

@dataclass
class PlayingCard:
    rank: str
    suit: str

@dataclass
class Deck:
    cards: List[PlayingCard]

>>> queen_of_hearts = PlayingCard('Q', 'Hearts')
>>> ace_of_spades = PlayingCard('A', 'Spades')
>>> two_cards = Deck([queen_of_hearts, ace_of_spades])


## Attempting to adjust example

@dataclass
class Protein:
    accession: str
    abundance_ratio: List[float]
    num_peptides: int
    peptides: List[peptides for Peptide if peptides.Protein == accession]


@dataclass
class Peptide:
    sequence: str
    abundance_ratio: List[float]



Values that can be added by calculations:
- av_abundance_ratio
- log10_pvalue
- gene_name
- gene_ontology

@dataclass
class Deck:
    cards: List[PlayingCard]
