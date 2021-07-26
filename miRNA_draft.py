import json


class Precursor():
    """A class with miRNA information
    """
    working_directory = "./"  # shared by all precursors

    def __init__(self, id_p, name, seq, org, ref, mirnas):  # variables unique to precursor
        self.precursor_ID = id_p
        self.precursor_name = name
        self.precursor_sequence = seq
        self.structure = ""  # z funkcji
        self.chromosome = []  # z funkcji
        self.genome_coordinates = []  # coordinantes of pre-miRNA start-end z funkcji
        self.strand = []  # z funkcji
        self.references = ref  # list of pubmed ids (RX line in miRNA.dat)
        self.organism = org  # full name
        self.taxonomy = []  # every tax level
        self.miRNAs = mirnas  # objects of miRNA class
        self.high_confidence = False

    def to_json(self):
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, indent=4)

    def dumper(self):
        try:
            return self.to_json()
        except:
            return self.__dict__

    def info(self):
        ref_join = '\n\t- '
        info = f"""
Precursor ID: {self.precursor_ID}
Precursor name: {self.precursor_name}
MiRNAs: {', '.join(self.miRNAs)}
Precursor structure: {self.structure}
Precursor sequence: {self.precursor_sequence}
Organism: {self.organism}
Taxonomy: {'/'.join(self.taxonomy)}
Chromosome: {self.chromosome}
Genome coordinates: {self.genome_coordinates}
Strand: {self.strand}
High confidence: {self.high_confidence}
References: {ref_join + ref_join.join(self.references)}
        """
        return info


class MiRNA():
    """A class with mature miRNA information
    """
    working_directory = "./"  # shared by all miRNAs

    def __init__(self, prec, m_id, name, pos, ev, ex, end):  # variables unique to miRNA
        self.precursor = prec  # może być ich kilka
        self.mature_name = name
        self.mature_ID = m_id
        self.pair_ID = None # ID of the other miRNA in pair
        self.mature_sequence = []
        self.mature_positions = [pos] # positions of mature miRNA on precursor in tab []
        self.evidence = ev
        self.experiment = ex
        self.end = end  # 3p, 5p, nothing
        self.chromosome_mi = []
        self.genome_coordinates_mi = []
        self.strand_mi = []
    
    def get_mature_seq(self, prec_seq, pos):
        seq = prec_seq[int(pos[0])-1:int(pos[1])]
        self.mature_sequence.append(seq)

    def to_json(self):
        return json.dumps(self, default=lambda o: o.__dict__, sort_keys=True, indent=4)

    def dumper(self):
        try:
            return self.to_json()
        except:
            return self.__dict__
