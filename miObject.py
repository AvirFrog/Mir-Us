"""
File with definition of objects used in the package. All of custom objects are defined here:
- Precursor object
- MiRNA object
"""

from collections import defaultdict
from colorama import init, Fore
import pprint as pp

__authors__ = ["Kacper Dudczak, Maciej Michalczyk"]
__copyright__ = "Copyright 2021, mirBase Project"
__credits__ = ["Marek Żywicki", "Marta Wysocka", "Kacper Dudczak", "Maciej Michalczyk"]
__license__ = "MIT"
__version__ = "0.2"
__maintainer__ = ["Kacper Dudczak", "Maciej Michalczyk"]
__email__ = ["kacper.dudczak19@gmail.com", "mccv99@gmail.com"]
__status__ = "Production"
__deprecated__ = False

# colorama setup
init(autoreset=True)


class Precursor:
    """
    Precursor class, which contains:
    - all data structures used to store information about precursors:
        - id
        - name
        - affiliated miRNAs
        - sequence
        - affiliated organism
        - full taxonomy of organism
        - genomic information:
            - chromosome
            - strand
            - coordinates
        - precursor confidence
        - references
    - utility functions:
        - pretty printing information from record (object)
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
        """
        Overrides print method to show Precursor object attributes in pretty and informative form.
        :return:  information about Precursor object attributes
        :rtype: str
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


#     def info(self):
#         """
#         Shows Precursor object attributes in pretty and informative form.
#         :return:  information about Precursor object attributes
#         :rtype: str
#         """
#
#         ref_join = '\n\t- '
#         info = f"""
# Precursor ID: {self.precursor_ID}
# Precursor name: {self.precursor_name}
# MiRNAs: {', '.join(self.miRNAs)}
# Precursor structure: {self.structure}
# Precursor sequence: {" " + self.precursor_sequence}
# Organism: {self.organism}
# Taxonomy: {'/'.join(self.taxonomy)}
# Chromosome: {''.join(self.chromosome)}
# Genome coordinates: {self.genome_coordinates}
# Strand: {''.join(self.strand)}
# High confidence: {self.high_confidence}
# References: {self.references}"""
#         return info
# References: {ref_join + ref_join.join(self.references)}


class MiRNA:
    """
    miRNA class, which contains:
    - all data structures used to store information about miRNA:
        - id
        - name
        - affiliated precursors
        - affiliated organism
        - mature sequence
        - genomic information:
            - chromosome
            - strand
            - coordinates
        - experiment backing
        - evidence
        - references
    - utility functions:
        - retrieving mature sequence
        - pretty printing information from record (object)
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
        """
        Overrides print method to show miRNA object attributes in pretty and informative form.
        :return:  information about miRNA object attributes
        :rtype: str
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
        """
        Returns miRNA's mature sequence
        :param prec_seq: precursor sequence
        :type prec_seq: str
        :param pos: positions of start and end for mature sequence
        :type pos: list
        :return: Returns full mature sequence of miRNA derived from precursor sequence
        :rtype: str
        """
        #return prec_seq[int(pos[0]) - 1:int(pos[1])]
        seq = prec_seq[int(pos[0]) - 1:int(pos[1])]
        self.mature_sequence.append(seq)

#     def info(self):
#         """
#         Shows miRNA object attributes in pretty and informative form.
#         :return:  information about miRNA object attributes
#         :rtype: str
#         """
#
#         ref_join = '\n\t- '
#         info = f"""
# Mature ID: {', '.join(self.mature_ID)}
# Mature name: {self.mature_name}
# Derivative: {self.precursor}
# Mature sequence: {self.mature_sequence}
# Mature positions: {self.mature_positions}
# Evidence: {self.evidence}
# Experiment: {self.experiment}
# End: {self.end}
# Chromosome: {self.chromosome_mi}
# Genome coordinates: {dict(self.genome_coordinates_mi)}
# Strand: {self.strand_mi}
# References: {self.references}"""
#         return info
