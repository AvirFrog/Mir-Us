"""
Provides definition of objects used in the package. All of custom objects are defined here:
    - Precursor object

    - MiRNA object
"""

from collections import defaultdict
from colorama import init, Fore
import pprint as pp

__authors__ = ["Kacper Dudczak, Maciej Michalczyk"]
__copyright__ = "Copyright 2021, mirBase Project"
__credits__ = ["Kacper Dudczak", "Maciej Michalczyk", "Marta Wysocka", "Marek Żywicki"]
__license__ = "MIT"
__version__ = "0.5"
__maintainer__ = ["Kacper Dudczak", "Maciej Michalczyk"]
__email__ = ["kacper.dudczak19@gmail.com", "mccv99@gmail.com"]
__status__ = "Production"
__deprecated__ = False

# colorama setup
init(autoreset=True)


class Precursor:
    """Precursor class to store information from a single Precursor record
    Attributes:
        precursor_ID (str): Precursor ID
        precursor_name (str): Precursor name
        miRNAs (list[str]): Affiliated miRNA IDs
        structure (str): Structure of precursor in dot-bracket format
        precursor_sequence (str): Nucleotide sequence of precursor
        organism (str): Name of affiliated organism
        taxonomy (list[str]): Full taxonomy of affiliated organism
        chromosome (list[str]): Chromosome names
        genome_coordinates (list[tuple(str, str)]): Genome coordinates of precursor
        strand (list[str]): Strand type ('+' or '-')
        high_confidence (bool): False by default
        references (list[str]): References from Pubmed (accession numbers)
    """

    def __init__(self, id, name, seq, org, ref, mirnas):
        self.precursor_ID = id
        self.precursor_name = name
        self.precursor_sequence = seq
        self.structure = ""  # dot-bracket
        self.chromosome = []
        self.genome_coordinates = []  # (start - end)
        self.strand = []  # '-' or '+'
        self.references = ref  # list of pubmed ids (RX line in miRNA.dat)
        self.organism = org  # full latin name
        self.taxonomy = []  # full taxonomy level
        self.miRNAs = mirnas  # objects of miRNA class
        self.high_confidence = False  # False by default

    def __repr__(self):
        """Overridden print method to show Precursor object attributes in pretty and informative form.

        Returns:
            str: All attributes from an object in 'Attribute: values' form.
        """
        info = f"""
        {Fore.YELLOW}Precursor ID: {Fore.RESET}{self.precursor_ID}
        {Fore.YELLOW}Precursor name: {Fore.RESET}{self.precursor_name}
        {Fore.YELLOW}MiRNAs: {Fore.RESET}{', '.join(self.miRNAs)}
        {Fore.YELLOW}Precursor structure: {Fore.RESET}{self.structure}
        {Fore.YELLOW}Precursor sequence: {Fore.RESET}{" " + self.precursor_sequence}
        {Fore.YELLOW}Organism: {Fore.RESET}{self.organism}
        {Fore.YELLOW}Taxonomy: {Fore.RESET}{'/'.join(self.taxonomy)}
        {Fore.YELLOW}Chromosome: {Fore.RESET}{', '.join(self.chromosome)}
        {Fore.YELLOW}Genome coordinates (start, end): {Fore.RESET}{', '.join(map(str, self.genome_coordinates))}
        {Fore.YELLOW}Strand: {Fore.RESET}{', '.join(self.strand)}
        {Fore.YELLOW}High confidence: {Fore.RESET}{self.high_confidence}
        {Fore.YELLOW}References: {Fore.RESET}{', '.join(map(str, self.references))}
"""
        return info


class MiRNA:
    """miRNA class to store information from a single miRNA record
    Attributes:
        mature_name (list[str]): miRNA names
        mature_ID (list[str]): miRNA IDs
        precursor (list[str]): IDs of affiliated precursors
        mature_sequence (list[str]): Mature miRNA sequences
        mature_positions (list[tuple(str, str)]): Pairs of mature seuqences positions from precursors
        organism (str): Name of affiliated organism
        evidence (list[str]): Evidence type
        experiment (list[str]): Types of conducted experiments to discover this miRNA (with numbers matching the reference index)
        end (list[str]): Strand end ('3p', '5p' or '-')
        chromosome_mi (list[str]): Chromosome names
        gen_coords (defaultdict[list[tuple(str, str)]]): Genome coordinates of miRNA from affiliated precursor
        strand_mi (list[str]): Strand type ('+' or '-')
        references (list[str]): References from Pubmed (accession numbers)
    """

    def __init__(self, prec, id, name, org, pos, evi, exp, end, ref):
        self.precursor = prec
        self.mature_name = name
        self.mature_ID = id
        self.organism = org
        self.pair_ID = None  # ID of the other miRNA in pair
        self.mature_sequence = []
        self.mature_positions = [pos]  # positions of mature miRNA on precursor in tab []
        self.evidence = evi
        self.experiment = exp
        self.end = end  # 3p, 5p, nothing
        self.chromosome_mi = []
        # self.genome_coordinates_mi = []
        self.genome_coordinates_mi = defaultdict(list)
        self.strand_mi = []
        self.references = ref

    def __repr__(self):
        """Overridden print method to show miRNA object attributes in pretty and informative form.

        Returns:
            str: All attributes from an object in 'Attribute: values' form.
        """
        gen_coords = pp.pformat(dict(self.genome_coordinates_mi))
        info = f"""
        {Fore.YELLOW}Mature ID: {Fore.RESET}{', '.join(self.mature_ID)}
        {Fore.YELLOW}Mature name: {Fore.RESET}{', '.join(self.mature_name)}
        {Fore.YELLOW}Derivative: {Fore.RESET}{', '.join(self.precursor)}
        {Fore.YELLOW}Mature sequence: {Fore.RESET}{', '.join(self.mature_sequence)}
        {Fore.YELLOW}Mature positions: {Fore.RESET}{', '.join(map(str, self.mature_positions))}
        {Fore.YELLOW}Organism: {Fore.RESET}{self.organism}
        {Fore.YELLOW}Evidence: {Fore.RESET}{', '.join(self.evidence)}
        {Fore.YELLOW}Experiment: {Fore.RESET}{', '.join(self.experiment)}
        {Fore.YELLOW}End: {Fore.RESET}{', '.join(self.end)}
        {Fore.YELLOW}Chromosome: {Fore.RESET}{', '.join(self.chromosome_mi)}
        {Fore.YELLOW}Genome coordinates: {Fore.RESET}{gen_coords}
        {Fore.YELLOW}Strand: {Fore.RESET}{', '.join(self.strand_mi)}
        {Fore.YELLOW}References: {Fore.RESET}{', '.join(map(str, self.references))}
"""
        return info

    def get_mature_seq(self, prec_seq, pos):
        """Sets miRNA's mature sequence

        Args:
            prec_seq (str): precursor sequence
            pos (list[str]): positions of start and end for mature sequence

        """
        # return prec_seq[int(pos[0]) - 1:int(pos[1])]
        seq = prec_seq[int(pos[0]) - 1:int(pos[1])]
        self.mature_sequence.append(seq)
