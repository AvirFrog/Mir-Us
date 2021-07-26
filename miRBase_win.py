import miRNA_draft
from collections import namedtuple as nt
from collections import defaultdict as dd
import argparse
import sys
import os
import gzip
import operator
from functools import reduce
import urllib.request
from functools import partial


class MiRBase():
    """ A class with information from miRBase
        and functions for this informations
    """
    def __init__(self, version):
        self.ftp_path = "ftp://mirbase.org/pub/mirbase/"
        self.miRBase_version = version
        self.working_directory = "./"
        self.miRNAs_ID = {}  # dct of miRNA ID = miRNA object
        self.precursors_ID = {}  # dct of pre-miRNA_ID = Precursor object
        self.org_sh = {}  # dct of 3-letter code and organism full name
        self.organisms = []  # list of tuples Organism = nt('Organism', 'organism division name tree taxid') a moze slownik?
        self.taxonomy_of_prec = dd(list)  # dct of every tax level with list of precursor and/or miRNAs IDs
        self.organisms_of_prec = dd(list) # dct of every organism with list of precursors and/or miRNAs IDs
        self.high_conf = []  #list of high confidance pre-miRNAs and/or miRNAs (IDs)
        self.done_coordinates = []  # list of 3-letter organisms codes when function get_coordinates for specific organism was performed
        self.structures = {}  # dct of pre-miRNA ID and their structures from miRNA.str file
        self.precursors_name = {}
        self.matures_name = {}
        self.system_name = ""
        self.pre_org = {}# key = 3 leter id of organism, value = list of prec
        self.mi_org = {}# key = 3 leter id of organism, value = list of mi

   

    def load_data(self):
        path = self.ftp_path + self.miRBase_version
        self.system_name = os.name

        name_organisms = path + "/organisms.txt.gz"
        self.load_organisms_txt(name_organisms)

        name_miRNAdat = path + "/miRNA.dat.gz"
        self.load_miRNAdat(name_miRNAdat)

        organism_codes = list(map(lambda x: getattr(x, "name"),self.organisms))  # lista wszystkich kodów organizmow
        tax_codes = list(map(lambda x: getattr(x, "tree"),self.organisms))
        tax_dct = dict(zip(organism_codes,tax_codes))

        self.get_hc()
        self.get_structures()

        for prec in self.precursors_ID:
            p_id = self.precursors_ID[prec].precursor_ID
            p_org = self.precursors_ID[prec].organism
            self.organisms_of_prec[p_org].append(p_id)
            tax_name = tax_dct[p_org].split(';')[:-1]
            self.precursors_ID[prec].taxonomy = tax_name
            for tax_level in tax_name:
                if not (p_id in self.taxonomy_of_prec[tax_level]):
                    self.taxonomy_of_prec[tax_level].append(p_id)
            self.precursors_ID[prec].structure = self.structures[prec]

    def load_organisms_txt(self, file_path):
        Organism = nt('Organism', 'organism division name tree taxid')
        with urllib.request.urlopen(file_path) as organisms_file_gz:
            with gzip.open(organisms_file_gz, mode = 'rt') as organisms_file:
                for line in organisms_file:
                    if not  line.startswith("#"):
                        tmp = line.split('\t')
                        self.org_sh[tmp[0]] = tmp[2]
                        org = Organism(tmp[0], tmp[1], tmp[2], tmp[3], tmp[4])
                        self.organisms.append(org)

    def load_miRNAdat(self, file_path):
        with urllib.request.urlopen(file_path) as miRNA_file_gz:
            with gzip.open(miRNA_file_gz, mode = 'rt') as miRNA_file:
                for line in miRNA_file:
                    if line.startswith("ID"):  
                        name_p = line.split()[1]
                        while not line.startswith("AC"):
                            line = next(miRNA_file)
                        id_p = line.split()[1][:-1]
                        while not line.startswith("DE"):
                            line = next(miRNA_file)
                        org = name_p.split("-")[0]
                        full_name = self.org_sh[org]
                        ref = []
                        while not line.startswith("FT"):
                            if line.startswith("RX"):
                                pubmed = line.split()[2][:-1]
                                ref.append(pubmed)
                            line = next(miRNA_file)
                        products = []
                        c = -1
                        while not line.startswith("SQ"):
                            next_line = True
                            if line.startswith("FT") and "miRNA" in line and ".." in line:
                                c += 1
                                start = line.split()[2].split("..")[0]
                                end = line.split()[2].split("..")[1].strip()
                                products.append([start,end])
                            elif line.startswith("FT") and "/accession=" in line:
                                ac = line.split("=")[1].replace('\"','').strip()
                                products[c].append(ac)
                            elif line.startswith("FT") and "/product=" in line:
                                name_mat = line.split("=")[1].replace('\"','').strip()
                                products[c].append(name_mat)
                            elif line.startswith("FT") and "/evidence=" in line:
                                ev = line.split("=")[1].replace('\"','').strip()
                                products[c].append(ev)
                                if ev == "experimental":
                                    line = next(miRNA_file)
                                    ex = ""
                                    if line.startswith("FT") and "/experiment=" in line:
                                        e = line.split("=")[1].replace('\"','').strip()
                                        ex += e
                                        if line.strip()[-1] != "\"":
                                            line = next(miRNA_file)
                                            next_line = False
                                            while line.startswith("FT") and "/experiment=" not in line and "miRNA" not in line:
                                                e = line.split("                   ")[1].replace('\"', '').strip()
                                                ex += " " + e
                                                line = next(miRNA_file)
                                        else:
                                            next_line = True
                                    products[c].append(ex)
                                else:
                                    products[c].append("-")
                            if next_line:
                                line = next(miRNA_file)
                        line = next(miRNA_file)
                        seq_full = ''
                        while not line.startswith("//"):
                            seq = line.split()[:-1]
                            seq_full += ''.join(seq)
                            line = next(miRNA_file)
                        m_ids = [m[2] for m in products]
                        if id_p not in self.precursors_ID:
                            prec = miRNA_draft.Precursor(id_p, name_p, seq_full, full_name, ref, m_ids)
                            self.precursors_ID[id_p] = prec
                            self.precursors_name[name_p] = id_p 
                        else:
                            for m in m_ids:
                                if m not in self.precursors_ID[id_p].miRNAs:
                                    self.precursors_ID[id_p].miRNAs.append(m)
                        for miRNA_vals in products:
                            if "5p" in miRNA_vals[3]:
                                e = "5p"
                            elif "3p" in miRNA_vals[3]:
                                e = "3p"
                            else:
                                e = "-"
                            if miRNA_vals[2] not in self.miRNAs_ID:
                                mature = miRNA_draft.MiRNA([id_p],[miRNA_vals[2]],[miRNA_vals[3]],[miRNA_vals[0],miRNA_vals[1]],[miRNA_vals[4]],[miRNA_vals[5]],[e])
                                mature.get_mature_seq(seq_full,[miRNA_vals[0],miRNA_vals[1]])
                                self.miRNAs_ID[miRNA_vals[2]] = mature
                            else:
                                self.miRNAs_ID[miRNA_vals[2]].precursor.append(id_p)
                                self.miRNAs_ID[miRNA_vals[2]].mature_name.append(miRNA_vals[3])
                                self.miRNAs_ID[miRNA_vals[2]].mature_positions.append([miRNA_vals[0],miRNA_vals[1]])
                                self.miRNAs_ID[miRNA_vals[2]].evidence.append(miRNA_vals[4])
                                self.miRNAs_ID[miRNA_vals[2]].experiment.append(miRNA_vals[5])
                                self.miRNAs_ID[miRNA_vals[2]].end.append(e)
                                self.miRNAs_ID[miRNA_vals[2]].get_mature_seq(seq_full, [miRNA_vals[0],miRNA_vals[1]])

    def get_coordinates(self, organism_code):
        if organism_code in self.done_coordinates:  # if not wystarczy !!!
            #get premi
            #return self.pre_org
            #get mi
            #return self.mi_org
            pass
        else:
            pre_org_list = []
            mi_org_list = []
            path = self.ftp_path + self.miRBase_version
            organism_c = path + '/genomes/' + organism_code + '.gff3'

            with urllib.request.urlopen(organism_c) as miRNA_file:
                file2 = miRNA_file.readlines()
            for line in file2:
                if not line.decode('UTF-8').startswith('#'):
                    split_name = line.decode('UTF-8').split('\t')
                    chr_miRNA = split_name[0]
                    miRNA_type = split_name[2]
                    start_seq = split_name[3].strip(' ')
                    end_seq = split_name[4].strip(' ')
                    strand_seq = split_name[6]
                    info_miRNA = split_name[8]
                    split_info_miRNA = info_miRNA.split(';')
                    sim_ID = split_info_miRNA[0].split('=')[1]
                    sim_Alias = split_info_miRNA[1].split('=')[1]
                    sim_name = split_info_miRNA[2].split('=')[1]

                    if miRNA_type == 'miRNA_primary_transcript':
                        # print(sim_Alias)
                        self.precursors_ID[sim_Alias].chromosome.append(chr_miRNA)
                        # print(self.precursors_ID[sim_Alias].chromosome)
                        self.precursors_ID[sim_Alias].strand.append(strand_seq)
                        # print(self.precursors_ID[sim_Alias].strand)
                        self.precursors_ID[sim_Alias].genome_coordinates.append([start_seq, end_seq])
                        # print(self.precursors_ID[sim_Alias].genome_coordinates)
                        # print('---' * 25)
                        # print(self.precursors_ID[sim_Alias].__dict__)
                        # print('---' * 25)
                        #self.pre_org[organism_code].append(sim_Alias)

                    elif miRNA_type == 'miRNA':
                        # !!!
                        # tutaj musisz wiedzieć na jakiej pozycji do tablicy to wpisać
                        # pozycja musi się zgadzać z pozycją pre-miRNA w zmiennej precursor
                        # do tego potrzebne będzie to pole Derives_from
                        # print(sim_Alias)
                        # self.miRNAs_ID[sim_Alias].chromosome_mi.append(chr_miRNA)
                        # print(self.miRNAs_ID[sim_Alias].chromosome_mi)
                        # self.miRNAs_ID[sim_Alias].strand_mi.append(strand_seq)
                        # print(self.miRNAs_ID[sim_Alias].strand_mi)
                        # self.miRNAs_ID[sim_Alias].genome_coordinates_mi.append([start_seq, end_seq])
                        # print(self.miRNAs_ID[sim_Alias].genome_coordinates_mi)
                        # # muszę znać indeks prekursora w tablicy prekursor
                        pass

        self.done_coordinates.append(organism_code)
        # print(self.pre_org)

        # uzupełnia obiekty miRNA o koordynaty genomowe (chr, start, stop) w oparciu o plik z adnotacjami dla konkretnego organizmu
        # 3 literowy skrotw nazwie pliku, rozpoznanie czy to miRNA czy prekursor
        # zabieram  start i stop nic i ID i chromosom
        #
        #
        #
        #
        #

    def get_structures(self):
        # na podstawie pliku .str tworzy uzupelnia obiekty Precursor o strukture w notacji dot-bracket
        iupac = ['a', 'c', 'g', 't', 'u', 'r', 'y', 's', 'w', 'k', 'm', 'b', 'd', 'h', 'v', 'n']
        start = 0
        end = 7
        x = 0
        path = self.ftp_path + self.miRBase_version
        name_str = path + "/miRNA.str.gz"
        
        with urllib.request.urlopen(name_str) as struct_file_gz:
            struct_file = gzip.open(struct_file_gz, mode = 'rt')
            #print("THE FILE IS BEING PROCESSED...")
            data = struct_file.readlines()
            lines_num = len(data)

            while x < lines_num:
                data_list = data[start:end]
                jdl = "".join(data_list).splitlines()
                new_list = [x for x in jdl if x]
                name = new_list[0]
                name2 = name.split()[0].strip('>')
                mrg_list = reduce(operator.add, zip(new_list[1], new_list[2]))
                line_3 = new_list[3].split()

                mrg_list2 = reduce(operator.add, zip(new_list[4], new_list[5]))

                # dot bracket
                all_dt_seq = ""
                join_mrg_list = "".join(mrg_list).replace(" ", "")

                dt_string1 = ""
                for letter, leter2 in zip(join_mrg_list, new_list[3]):
                    if (letter.lower() in iupac) and (leter2 == '|'):
                        dt_string1 = dt_string1 + '('
                    elif (letter.lower() in iupac) and not (leter2 == '|'):
                        dt_string1 = dt_string1 + '.'
                    else:
                        continue

                acgu_dt = [x.lower() for x in line_3 if x.lower() in iupac]
                if len(acgu_dt) == 1:
                    acgu_dt = ''.join(acgu_dt).lower().replace("-", "")#.replace('a', '.').replace('c', '.').replace('g','.').replace(
                        #'u', '.').replace('n', '.').replace('t', '.').replace('r', '.').replace('y','.').replace('s', '.').replace(
                        #'w', '.').replace('k', '.').replace('m', '.').replace('b', '.').replace('d', '.').replace('h', '.').replace(
                        #'v', '.')
                    for nt in iupac:
                        if nt in acgu_dt:
                            acgu_dt = acgu_dt.replace(nt,'.')
                    dt_string1 = dt_string1 + acgu_dt

                join_mrg_list2 = "".join(mrg_list2).replace(" ", "")
                dt_string2 = ""

                for letter, leter2 in zip(join_mrg_list2, new_list[3]):
                    if (letter.lower() in iupac) and (leter2 == '|'):
                        dt_string2 = dt_string2 + ')'
                    elif (letter.lower() in iupac) and not (leter2 == '|'):
                        dt_string2 = dt_string2 + '.'

                revers_dt_string2 = dt_string2[::-1]

                all_dt_seq = dt_string1 + revers_dt_string2

                start = start + 8
                end = end + 8
                x = x + 8

                id_prec = self.precursors_name[name2]
                self.structures[id_prec] = all_dt_seq

            struct_file.close()
        #print("FILE PROCESSING COMPLETED.")

    def get_hc(self):
        path = self.ftp_path + self.miRBase_version
        name_hc = path + "/hairpin_high_conf.fa.gz"
        
        with urllib.request.urlopen(name_hc) as hc_file_gz:
            with gzip.open(hc_file_gz, mode = 'rt') as hc_file:
                for line1 in hc_file:
                    if line1.startswith('>'):
                        tmp_hc =line1.split()[1]
                        self.precursors_ID[tmp_hc].high_confidance = True
                        self.high_conf.append(tmp_hc)

    def get_mature_miRNA(self, id = "", name = "", organism_name = "", organism_code = "", tax_level = "", chr = "", start = "", end = "", precursor_ID = ""):
        # zwraca liste nazw miRNA dla danego organizmu, poziomu taksonomicznego
        # zwraca nazwe miRNA dla danego ID
        # zwraca ID miRNA dla danej nazwy
        # zwaraca liste ID miRNA dla podanych koordynatów (organizm, chromosom, start, stop - wszystkie 4 musza być podane) - do wykorzystania funckja get_coordinates
        # zwaca listę ID miRNA dla danego ID prekursora
        pass

    def get_precursor(self, id="", name="", organism_name="", organism_code="", tax_level="", chr="", start="",
                      end="", strand='', mirna_ID=""):

        if organism_code:
            self.get_coordinates(organism_code)
        else:
            for key, value in self.org_sh.items():
                if organism_name == value:
                    self.get_coordinates(key)

        # !!! funkcja get_coordinates ma byc wywolana jesli ktos poda
        # 1) (organism_name lub organism_code) oraz chr
        # 2) (organism_name lub organism_code) oraz chr, start i end
        # 3) (organism_name lub organism_code) oraz strand
        # 4) (organism_name lub organism_code) oraz chr i strand
        # 5) (organism_name lub organism_code) oraz chr, strand, start i end

        if id != '':
            # !!! przeciez self.precursors_ID[id].precursor_name zawiera nazwe, nie ma tu co szukac forem, wystarczy to zwrocic
            try:
                name_p = self.precursors_ID[id].precursor_name
                print(name_p)
            except KeyError:
                print('ErRoR')

        elif name != '':
            id_pr = self.precursors_name[name]
            return id_pr
            # print(self.precursors_name[name])

        elif (organism_code != '' or organism_name != '') and (chr == '') and (start == '' and end == ''):
            # !!! masz slownik gdzie dla kazdego organizmu sa podane ID pre-miRNA ...
            # !!! zeby otrzymac klucz dla danej wartosci ze slownika mozna zrobic tak:
            # !!! code = list(self.org_sh.keys())[list(self.org_sh.values()).index(organism_name)]
            if organism_code:
                id_pr_list = []
                for key in self.precursors_name.keys():
                    if organism_code in key:
                        id_pr = self.precursors_name[key]
                        id_pr_list.append(id_pr)
                return id_pr_list
                # print(len(id_pr_list))
            else:
                # self.org_sh

                for key, value in self.org_sh.items():
                    if organism_name == value:
                        code = key
                        id_pr_list = []
                        for key2 in self.precursors_name.keys():
                            if code in key2:
                                id_pr = self.precursors_name[key2]
                                id_pr_list.append(id_pr)
                        return id_pr_list
                        # print(id_pr_list)
                        # print(len(id_pr_list))

        elif (organism_code != '' or organism_name != '') and (chr != '') and (start == '' and end == ''):
            if organism_code:
                naame = self.org_sh[organism_code]
                id_pr_list = []
                for value in self.precursors_ID.values():
                    if chr in value.chromosome:
                        if naame == value.organism:
                            id_pr_list.append(value.precursor_ID)
                return id_pr_list
                # print(id_pr_list)
                # print(len(id_pr_list))
            else:
                id_pr_list = []
                for value in self.precursors_ID.values():
                    if chr in value.chromosome:
                        if organism_name == value.organism:
                            id_pr_list.append(value.precursor_ID)
                return id_pr_list
                # print(id_pr_list)
                # print(len(id_pr_list))


        elif (organism_code != '' or organism_name != '') and (chr != '') and (start != '' and end != ''):

            if organism_code:
                naame = self.org_sh[organism_code]
                id_pr_list = []
                for value in self.precursors_ID.values():
                    if chr in value.chromosome:
                        if naame == value.organism:
                            for x in value.genome_coordinates:
                                if (int(start) <= int(x[0])) and (int(end) >= int(x[1])):
                                    # !!! start >= x[0] end <= x[1]
                                    id_pr_list.append(value.precursor_ID)
                return id_pr_list
                # print(id_pr_list)
                # print(len(id_pr_list))
            else:
                id_pr_list = []
                for value in self.precursors_ID.values():
                    if chr in value.chromosome:
                        if organism_name == value.organism:
                            for x in value.genome_coordinates:
                                if (int(start) <= int(x[0])) and (int(end) >= int(x[1])):
                                    id_pr_list.append(value.precursor_ID)
                return id_pr_list
                # print(id_pr_list)
                # print(len(id_pr_list))

        elif mirna_ID != '':
            pr_id = self.miRNAs_ID[mirna_ID].precursor
            # print(pr_id)
            return pr_id
        elif tax_level != '':
            list_of_premi = self.taxonomy_of_prec[tax_level]
            return list_of_premi
            # print(self.taxonomy_of_prec[tax_level])

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--version", help = "database version", required = True, type = str, dest ="v")
    args = parser.parse_args()
    mirbase = MiRBase(version=args.v)
    mirbase.load_data()
    #mirbase.get_coordinates(organism_code='aga')
 #ftp://mirbase.org/pub/mirbase/CURRENT/genomes/aga.gff3
    mirbase.get_precursor(id='MI0000001')