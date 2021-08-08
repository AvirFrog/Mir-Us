"""
Main file of the package. All of the actions are handled with this code:
- parsing mirBase files and compiling them to .mir files
- reading data from .mir files
- accessing data from .mir file by package user
"""

import functools
import gzip
import json
import operator
import urllib.request
from collections import defaultdict as dd
from collections import namedtuple as nt
from functools import reduce
from timeit import default_timer as timer
from varname import nameof
import numpy as np
import pprint as pp

import dill
from colorama import init, Fore

import miObject

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


# DECORATORS
def time_this(func):
    """
    Decorator which returns information about execution of decorated function
    :param func: Any miBase function
    :return: Execution time and values returned by a function
    """

    @functools.wraps(func)
    def wrapper_timer(*args, **kwargs):
        start = timer()
        values = func(*args, **kwargs)
        end = timer()
        runtime = end - start
        if values is None:
            print(f"{Fore.RED}[Mir-Us]  {func.__name__!r} No records matching given criteria.")
        else:
            print(f"{Fore.GREEN}[Mir-Us]  {func.__name__!r} found {values[1]} results in {runtime:.6f} seconds")
            return values[0]

    return wrapper_timer


# CLASSES
class MiRBase:
    """
    Main class, which contains:
    - all data structures used to store information
    - all functions to retrieve information by user
    - all functions to load data from indexed files (.mir files)
    - all functions to compile data into indexed files (.mir files)
    """

    def __init__(self, version="CURRENT"):
        self._ftp_path = "ftp://mirbase.org/pub/mirbase/"
        self._miRBase_version = version
        self._working_directory = "./"
        self._miRNAs_ID = {}  # dct of miRNA ID = miRNA object
        self._precursors_ID = {}  # dct of pre-miRNA_ID = Precursor object
        self._org_sh = {}  # dct of 3-letter code and organism full name
        self._organisms = []  # list of tuples Organism = nt('Organism', 'organism division name tree taxid') a moze slownik?
        self._taxonomy_of_prec = dd(list)  # dct of every tax level with list of precursor and/or miRNAs IDs
        self._organisms_of_prec = dd(list)  # dct of every organism with list of precursors and/or miRNAs IDs
        self._high_conf = []  # list of high confidence pre-miRNAs and/or miRNAs (IDs)
        self._done_coordinates = []  # list of 3-letter organisms codes when function get_coordinates for specific organism was performed
        self._structures = {}  # dct of pre-miRNA ID and their structures from miRNA.str file
        self._precursors_name = {}
        self._matures_name = {}
        self._system_name = ""
        self._pre_org = {}  # key = 3 leter id of organism, value = list of prec
        self._mi_org = {}  # key = 3 leter id of organism, value = list of mi
        self._prec_genome = {}
        self._miRNA_genome = {}

        self._loader = MiRLoad  # reference to loader class - very important

        self._Organism = nt('Organism', 'organism division name tree taxid')

        try:
            self._load_all_data()
        except:
            print(f"{Fore.RED}[Mir-Us]   Missing data files (.mir); "
                  f"{Fore.YELLOW}performing compilation of data files...")
            # try:
            self._compile_indexes()
            print(f"{Fore.YELLOW}[Mir-Us]   Data files compiled successfully!")
            self._load_all_data()
            # except:
            #     print(f"{Fore.RED}[Mir-Us]   Error! Cannot compile data files.")
            #     return

    def _compile_indexes(self):
        """
        Makes .mir files which contain indexed data.
        :return: .mir files
        """
        path = self._ftp_path + self._miRBase_version

        # COMPILE ORGANISMS--------------------------------------
        self._loader.load_organisms(self, file_path=path + "/organisms.txt.gz")
        with open('organisms.mir', 'wb') as fh_org_dump:
            dill.dump(self._organisms, fh_org_dump)
        fh_org_dump.close()

        with open('org_short.mir', 'wb') as fh_orgsh_dump:
            dill.dump(self._org_sh, fh_orgsh_dump)
        fh_orgsh_dump.close()
        # -------------------------------------------------------

        # COMPILE MIRNA------------------------------------------
        self._loader.load_miRNA(self, file_path=path + "/miRNA.dat.gz")
        self._loader.load_genome(self, file_path=path + "/genomes/")
        with open('precursors_ID.mir', 'wb') as fh_precid_dump:
            dill.dump(self._precursors_ID, fh_precid_dump)
        fh_precid_dump.close()

        with open('precursors_name.mir', 'wb') as fh_precname_dump:
            dill.dump(self._precursors_name, fh_precname_dump)
        fh_precname_dump.close()

        with open('miRNAs_ID.mir', 'wb') as fh_mirnasid_dump:
            dill.dump(self._miRNAs_ID, fh_mirnasid_dump)
        fh_mirnasid_dump.close()

        # with open('org_short.mir', 'wb') as fh_orgsh_dump:
        #     dill.dump(self._org_sh, fh_orgsh_dump)
        # fh_orgsh_dump.close()
        # -------------------------------------------------------

        # COMPILE HIGH-CONF--------------------------------------
        self._loader.load_hc(self, file_path=path + "/hairpin_high_conf.fa.gz")
        with open('high_conf.mir', 'wb') as fh_high_dump:
            dill.dump(self._high_conf, fh_high_dump)
        fh_high_dump.close()
        # -------------------------------------------------------

        # COMPILE STRUCTURES-------------------------------------
        self._loader.load_structures(self, file_path=path + "/miRNA.str.gz")
        with open('structures.mir', 'wb') as fh_mirstruc_dump:
            dill.dump(self._structures, fh_mirstruc_dump)
        fh_mirstruc_dump.close()
        # -------------------------------------------------------

        # COMPILE TAXONOMY---------------------------------------
        self._loader.load_taxonomy(self)
        with open('taxonomy_prec.mir', 'wb') as fh_taxprec_dump:
            dill.dump(self._taxonomy_of_prec, fh_taxprec_dump)
        fh_taxprec_dump.close()

        with open('taxonomy_org.mir', 'wb') as fh_taxorg_dump:
            dill.dump(self._organisms_of_prec, fh_taxorg_dump)
        fh_taxorg_dump.close()
        # -------------------------------------------------------

    def _load_all_data(self):
        """
        Loads all data from .mir files
        :return: All data structures have their data assigned.
        """
        # LOAD ORGANISMS-----------------------------------------
        with open('organisms.mir', 'rb') as fh_org_load:
            self._organisms = dill.load(fh_org_load)
        fh_org_load.close()
        # -------------------------------------------------------

        # LOAD MIRNA---------------------------------------------
        with open('precursors_ID.mir', 'rb') as fh_precid_load:
            self._precursors_ID = dill.load(fh_precid_load)
        fh_precid_load.close()

        with open('precursors_name.mir', 'rb') as fh_precname_load:
            self._precursors_name = dill.load(fh_precname_load)
        fh_precid_load.close()

        with open('miRNAs_ID.mir', 'rb') as fh_mirnasid_load:
            self._miRNAs_ID = dill.load(fh_mirnasid_load)
        fh_mirnasid_load.close()

        with open('org_short.mir', 'rb') as fh_orgsh_load:
            self._org_sh = dill.load(fh_orgsh_load)
        fh_orgsh_load.close()
        # -------------------------------------------------------

        # LOAD HIGH-CONF-----------------------------------------
        with open('high_conf.mir', 'rb') as fh_high_load:
            self._high_conf = dill.load(fh_high_load)
        fh_high_load.close()
        # -------------------------------------------------------

        # LOAD STRUCTURES----------------------------------------
        with open('structures.mir', 'rb') as fh_mirstruc_load:
            self._structures = dill.load(fh_mirstruc_load)
        fh_mirstruc_load.close()
        # -------------------------------------------------------

        # LOAD TAXONOMY------------------------------------------
        with open('taxonomy_prec.mir', 'rb') as fh_taxprec_load:
            self._taxonomy_of_prec = dill.load(fh_taxprec_load)
        fh_taxprec_load.close()
        with open('taxonomy_org.mir', 'rb') as fh_taxorg_load:
            self._organisms_of_prec = dill.load(fh_taxorg_load)
        fh_taxorg_load.close()
        # -------------------------------------------------------

        # PERFORM MERGE------------------------------------------
        self._merge_data()
        # -------------------------------------------------------

    def _merge_data(self):
        """
        Merges data - some data structures lacks certain parts and this function completes it.
        :return: Data structures completed with information which is connected between them.
        """
        # MERGE HIGH_CONF----------------------------------------
        for item in self._high_conf:
            self._precursors_ID[item].high_confidence = True
        # -------------------------------------------------------

        # MERGE STRUCTURES---------------------------------------
        for prec in self._precursors_ID:
            self._precursors_ID[prec].structure = self._structures[prec]
        # -------------------------------------------------------

        # MERGE TAXONOMY-----------------------------------------
        organism_codes = list(map(lambda x: getattr(x, "name"), self._organisms))  # list of all organism codes
        tax_codes = list(map(lambda x: getattr(x, "tree"), self._organisms))
        tax_codes = [x.rstrip(";").split(";") for x in tax_codes]
        tax_dct = dict(zip(organism_codes, tax_codes))

        for prec in self._precursors_ID:
            self._precursors_ID[prec].taxonomy = tax_dct[self._precursors_ID[prec].organism]
        # -------------------------------------------------------

    def _precursor_retrieve(self, passed_id, current_result):
        try:
            if not self._exists(current_result, self._precursors_ID[passed_id]):
                return self._precursors_ID[passed_id]
        except KeyError as e:
            pass

    def _mirna_retrieve(self, passed_id, current_result):
        try:
            if not self._exists(current_result, self._miRNAs_ID[passed_id]):
                return self._miRNAs_ID[passed_id]
        except KeyError as e:
            pass

    def _input(self, prompt, type=int):
        while True:
            try:
                return type(input(prompt))
            except:
                pass
        return res

    def get_organisms_list(self):
        """
        Returns list of organisms.
        :return: list of organism namedtuples
        :rtype: list
        """
        return self._organisms

    @time_this
    def get_tax_level(self, organism=None):
        """
        Returns taxonomy level assigned to organism.
        :param organism: list of names (strings) of organisms
        :type organism: list
        :return: dictionary of organisms and its assigned taxonomy (key: organism, value: taxonomy)
        :rtype: dict
        """
        if isinstance(organism, str):
            organism = [organism]
        result = {}
        organism_codes = list(map(lambda x: getattr(x, "name"), self._organisms))  # lista wszystkich kodów organizmow
        tax_codes = list(map(lambda x: getattr(x, "tree"), self._organisms))
        tax_codes = [x.rstrip(";").split(";") for x in tax_codes]
        tax_dct = dict(zip(organism_codes, tax_codes))
        try:
            for org in organism:
                try:
                    result[org] = tax_dct[org]
                except:
                    continue
        except:
            pass
        if not result:
            return None
        return result, len(result)

    @time_this
    def get_organism(self, tax=None):
        """
        Returns organisms which are assigned to a given taxonomy level.
        :param tax: string representing taxonomy level
        :type tax: str
        :return: list of strings representing organisms
        :rtype: list
        """
        result = []
        organism_codes = list(map(lambda x: getattr(x, "name"), self._organisms))  # lista wszystkich kodów organizmow
        tax_codes = list(map(lambda x: getattr(x, "tree"), self._organisms))
        tax_codes = [x.rstrip(";").split(";") for x in tax_codes]
        tax_dct = dict(zip(organism_codes, tax_codes))
        #print(tax_dct)
        for key, value in tax_dct.items():
            if tax in value:
                result.append(key)
            else:
                continue
        if not result:
            return None
        return result, len(result)

    def get_organisms_short(self):
        """
        Returns dictionary of all organism name and its abbreviation
        :return: dictionary (key: abbreviation, value: full organism name)
        :rtype: dict
        """
        return self._org_sh

    @time_this
    def get_taxid(self, organism=None):
        """
        Returns taxid assigned to organism
        :param organism: list of strings representing organisms
        :type organism: list
        :return: dictionary of organisms with assigned taxid (key: full organism name, value: taxid)
        :rtype: dict
        """
        if isinstance(organism, str):
            organism = [organism]
        result = {}
        organism_codes = list(map(lambda x: getattr(x, "name"), self._organisms))  # lista wszystkich kodów organizmow
        tax_id = list(map(lambda x: getattr(x, "taxid"), self._organisms))
        tax_dct = dict(zip(organism_codes, tax_id))
        try:
            for org in organism:
                try:
                    result[org] = tax_dct[org]
                except:
                    continue
        except:
            pass
        if not result:
            return None
        return result, len(result)

    @time_this
    def get_precursor(self, id: list = None, name="", organism_name="", tax_level="", chr="", start="",
                      end="", strand='', mirna_id=""):
        """
        Returns precursor objects according to given search criteria
        :param id: list of strings representing ids of precursors
            (default is None)
        :type id: list
        :param name: name of precursor
        :type name: str
        :param organism_name: organism name in which precursor is present
        :type organism_name: str
        :param tax_level: taxonomy level affiliated with precursor
        :type tax_level: str
        :param chr: chromosome name in which precursor is present
        :type chr: str
        :param start: genomic location in which precursor sequence starts
        :type start: str
        :param end: genomic location in which precursor sequence ends
        :type end: str
        :param strand: strand name in which precursor is present
        :type strand: str
        :param mirna_id: string representing id of miRNA affiliated with precursor
        :type mirna_id: str
        :return: list of Precursor objects
        :rtype: list
        """
        if isinstance(id, str):
            id = [id]
        result = []
        if id:
            result = [self._precursor_retrieve(i, result) for i in id if
                      self._precursor_retrieve(i, result) is not None]
        if name:
            for prec in self._precursors_ID:
                if (self._precursors_ID[prec].precursor_name == name) and not self._exists(result,
                                                                                            self._precursors_ID[prec]):
                    result.append(self._precursors_ID[prec])
        if organism_name:
            for prec in self._precursors_ID:
                if (self._precursors_ID[prec].organism == organism_name) and not self._exists(result,
                                                                                               self._precursors_ID[
                                                                                                   prec]):
                    result.append(self._precursors_ID[prec])
        if tax_level:
            for prec in self._precursors_ID:
                if (tax_level in self._precursors_ID[prec].taxonomy) and not self._exists(result,
                                                                                           self._precursors_ID[prec]):
                    result.append(self._precursors_ID[prec])
        if mirna_id:
            for elem in mirna_id:
                try:
                    for prec in self._miRNAs_ID[elem].precursor:
                        if not self._exists(result, prec):
                            result.append(self._precursors_ID[prec])
                except:
                    continue
        if chr:
            for prec in self._precursors_ID:
                if (chr in self._precursors_ID[prec].chromosome) and not self._exists(result,
                                                                                       self._precursors_ID[prec]):
                    result.append(self._precursors_ID[prec])
        if strand:
            for prec in self._precursors_ID:
                if (strand in self._precursors_ID[prec].strand) and not self._exists(result,
                                                                                      self._precursors_ID[prec]):
                    result.append(self._precursors_ID[prec])
        if start and end:
            int_start = int(start)
            int_end = int(end)
            if not int_start < int_end:
                print("{Fore.RED}[Mir-Us]   Wrong coordinates; start cannot be lower than end")
                return
            for prec in self._precursors_ID:
                try:
                    f_start = int(self._precursors_ID[prec].genome_coordinates[0][0])
                    f_stop = int(self._precursors_ID[prec].genome_coordinates[0][1])
                except:
                    continue
                if (int_start <= f_start < int_end) and (int_start < f_stop <= int_end):
                    result.append(self._precursors_ID[prec])
        if not result:
            return None

        return result, len(result)

    @time_this
    def get_references(self, mirna_id=None, mirna_name=None, prec_id=None, link=False):
        """
        Returns list of reference numbers
        :param prec_id: list of strings representing ids of precursors
        :param link: bool value which forces function to return PubMed links instead of reference numbers
        :param mirna_id: list of strings representing ids of miRNA
            (default is None)
        :type mirna_id: list
        :param mirna_name: list of strings representing names of miRNA
            (default is None)
        :type mirna_name: list
        :return: list of strings representing reference numbers
        :rtype: list
        """
        if isinstance(mirna_id, str):
            mirna_id = [mirna_id]
        if isinstance(mirna_name, str):
            mirna_name = [mirna_name]
        if isinstance(prec_id, str):
            prec_id = [prec_id]
        result = []
        if mirna_id:
            for m_id in mirna_id:
                try:
                    for ref in self._miRNAs_ID[m_id].references:
                        if ref not in result and not link:
                            result.append(ref)
                        elif link:
                            result.append(f"https://pubmed.ncbi.nlm.nih.gov/{ref}/")
                except:
                    continue
        if mirna_name:
            for mi_name in mirna_name:
                try:
                    for mi in self._miRNAs_ID:
                        if mi_name in self._miRNAs_ID[mi].mature_name:
                            for ref in self._miRNAs_ID[mi].references:
                                if ref not in result and not link:
                                    result.append(ref)
                                elif link:
                                    result.append(f"https://pubmed.ncbi.nlm.nih.gov/{ref}/")

                except:
                    continue
        if prec_id:
            for p_id in prec_id:
                try:
                    for ref in self._precursors_ID[p_id].references:
                        if ref not in result and not link:
                            result.append(ref)
                        elif link:
                            result.append(f"https://pubmed.ncbi.nlm.nih.gov/{ref}/")
                except:
                    continue
        if not result:
            return None
        return result, len(result)

    @time_this
    def get_structure(self, id=None, name=None):
        """
        Returns dictionary of precursors ids with assigned structures in dot-bracket format
        :param name: list of strings representing precursors name
            (default is None)
        :param id: list of strings representing precursors ids
            (default is None)
        :type id: list
        :return: dictionary of precursors ids with assigned structure (key: precursor id, value: structure)
        :rtype: dict
        """
        result = {}
        if isinstance(id, str):
            id = [id]
        if isinstance(name, str):
            name = [name]
        if id:
            try:
                for i in id:
                    try:
                        result[self._precursors_ID[i].precursor_ID] = self._precursors_ID[i].structure
                    except:
                        continue
            except:
                pass
        if name:
            try:
                for i in name:
                    try:
                        result = {self._precursors_ID[prec].precursor_name: self._precursors_ID[prec].structure for
                                  prec in self._precursors_ID if self._precursors_ID[prec].precursor_name == i}
                    except:
                        continue
            except:
                pass
        if not result:
            return None
        return result, len(result)

    def _exists(self, to_compare, obj):
        """
        Utility function, which compares list of objects to particular object, to check if compared object is present
        in the list
        :param to_compare: list of objects of any type
        :type to_compare: list
        :param obj: object of any type to be compared with list of objects
        :type obj: object
        :return: returns True if compared object is present in the list of objects
        :rtype: bool
        """
        for elem in to_compare:
            if elem is obj:
                return True

    @time_this
    def get_mirna(self, mirna_id: list = None, name="", organism_name="", tax_level="", chr="", start="",
                  end="", strand='', prec_id: list = None):
        """
        Returns miRNA objects according to given search criteria
        :param prec_id: list of strings representing ids of precursors
            (default is None)
        :type prec_id: list
        :param mirna_id: list of strings representing ids of miRNAs
            (default is None)
        :type mirna_id: list
        :param name: name of miRNA
        :type name: str
        :param organism_name: organism name in which miRNA is present
        :type organism_name: str
        :param tax_level: taxonomy level affiliated with miRNA
        :type tax_level: str
        :param chr: chromosome name in which miRNA is present
        :type chr: str
        :param start: genomic location in which miRNA sequence starts
        :type start: str
        :param end: genomic location in which miRNA sequence ends
        :type end: str
        :param strand: strand name in which miRNA is present
        :type strand: str
        :return: list of miRNA objects
        :rtype: list
        """
        if isinstance(mirna_id, str):
            mirna_id = [mirna_id]
        #print(locals())
        dict_result = {}
        #result = []
        first = True

        if not end and start or start and end or organism_name and chr or organism_name and strand:
            result = []
            if organism_name:
                for mi in self._miRNAs_ID:
                    if (self._miRNAs_ID[mi].organism == organism_name) and not self._exists(result, self._miRNAs_ID[mi]):
                        result.append(self._miRNAs_ID[mi])
                first = False
            if chr:
                temp_result = []
                if first:
                    for mi in self._miRNAs_ID:
                        if (chr in self._miRNAs_ID[mi].chromosome_mi) and not self._exists(result, self._miRNAs_ID[mi]):
                            result.append(self._miRNAs_ID[mi])
                    first = False
                elif not first:
                    #count = 0
                    #print(len(result))
                    for res in result:
                        #count += 1
                        #print(f"Current {(count)}: {res.chromosome_mi}")
                        if chr in res.chromosome_mi and not self._exists(temp_result, res):
                            #print(f"Kept:{res.chromosome_mi}")
                            temp_result.append(res)
                result = temp_result
                #print("Chr from genomic context")
            if strand:
                temp_result = []
                if first:
                    for mi in self._miRNAs_ID:
                        if (strand in self._miRNAs_ID[mi].strand_mi) and not self._exists(result, self._miRNAs_ID[mi]):
                            result.append(self._miRNAs_ID[mi])
                    first = False
                elif not first:
                    for res in result:
                        if strand in res.strand_mi and not self._exists(temp_result, res):
                            temp_result.append(res)
                result = temp_result
                #print("Strand from genomic context")
            if start and end:
                temp_result = []
                if first:
                    int_start = int(start)
                    int_end = int(end)
                    if not int_start < int_end:
                        print(f"{Fore.RED}[Mir-Us]   Wrong coordinates; start cannot be lower than end")
                        return None
                    for mi in self._miRNAs_ID:
                        for key in self._miRNAs_ID[mi].genome_coordinates_mi:
                            for coord in self._miRNAs_ID[mi].genome_coordinates_mi[key]:
                                try:
                                    f_start = int(coord[0])
                                    f_stop = int(coord[1])
                                except:
                                    continue
                                if (int_start <= f_start < int_end) and (int_start < f_stop <= int_end)\
                                        and not self._exists(result, self._miRNAs_ID[mi]):
                                    result.append(self._miRNAs_ID[mi])
                    first = False
                elif not first:
                    int_start = int(start)
                    int_end = int(end)
                    if not int_start < int_end:
                        print(f"{Fore.RED}[Mir-Us]   Wrong coordinates; start cannot be lower than end")
                        return None
                    for res in result:
                        for key in res.genome_coordinates_mi:
                            for coord in res.genome_coordinates_mi[key]:
                                try:
                                    f_start = int(coord[0])
                                    f_stop = int(coord[1])
                                except:
                                    continue
                                if (int_start <= f_start < int_end) and (int_start < f_stop <= int_end)\
                                        and not self._exists(temp_result, res):
                                    temp_result.append(res)
                    result = temp_result
            if not end and start:
                temp_result = []
                if first:
                    int_start = int(start)
                    for mi in self._miRNAs_ID:
                        for key in self._miRNAs_ID[mi].genome_coordinates_mi:
                            for coord in self._miRNAs_ID[mi].genome_coordinates_mi[key]:
                                try:
                                    f_start = int(coord[0])
                                except:
                                    continue
                                if (int_start <= f_start) and not self._exists(
                                        result, self._miRNAs_ID[mi]):
                                    result.append(self._miRNAs_ID[mi])
                    first = False
                elif not first:
                    int_start = int(start)
                    #count = 0
                    for res in result:
                        #count += 1
                        for key in res.genome_coordinates_mi:
                            for coord in res.genome_coordinates_mi[key]:
                                try:
                                    f_start = int(coord[0])
                                    #print(f"Current {(count)}: {f_start}")
                                except:
                                    continue
                                if int_start <= f_start and not self._exists(temp_result, res):
                                    temp_result.append(res)
                                    #print(f"Kept:{f_start}")
                    result = temp_result
                #print("Only start")
            dict_result["genomic-search"] = result
        if mirna_id:
            result = []
            result = [self._mirna_retrieve(i, result) for i in mirna_id if self._mirna_retrieve(i, result) is not None]
            dict_result["mirna_id-search"] = result
        if prec_id:
            result = []
            for elem in prec_id:
                    try:
                        for mi in self._precursors_ID[elem].miRNAs:
                            if not self._exists(result, mi):
                                result.append(self._miRNAs_ID[mi])
                    except:
                        continue
            dict_result["prec_id-search"] = result
        if name:
            result = []
            for mi in self._miRNAs_ID:
                if (self._miRNAs_ID[mi].mature_name == name) and not self._exists(result, self._miRNAs_ID[mi]):
                    result.append(self._miRNAs_ID[mi])
            dict_result["name-search"] = result
        if not start and not end and not chr and not strand and organism_name:
            result = []
            for mi in self._miRNAs_ID:
                if (self._miRNAs_ID[mi].organism == organism_name) and not self._exists(result, self._miRNAs_ID[mi]):
                    result.append(self._miRNAs_ID[mi])
            dict_result["organism-search"] = result
            #print("Just org")
        if tax_level:
            result = []
            for prec in self._precursors_ID:
                if tax_level in self._precursors_ID[prec].taxonomy:
                    mi_id = self._precursors_ID[prec].miRNAs
                    for elem in mi_id:
                        if not self._exists(result, self._miRNAs_ID[elem]):
                            result.append(self._miRNAs_ID[elem])
            dict_result["tax-search"] = result
        if not start and not end and not organism_name and not strand and chr:
            result = []
            for mi in self._miRNAs_ID:
                if (chr in self._miRNAs_ID[mi].chromosome_mi) and not self._exists(result, self._miRNAs_ID[mi]):
                    result.append(self._miRNAs_ID[mi])
            dict_result["chr-search"] = result
            #print("Chr by mistake but will do the job either way")
        if not start and not end and not organism_name and not chr and strand:
            result = []
            for mi in self._miRNAs_ID:
                if (strand in self._miRNAs_ID[mi].strand_mi) and not self._exists(result, self._miRNAs_ID[mi]):
                    result.append(self._miRNAs_ID[mi])
            dict_result["strand-search"] = result
            #print("Strand by mistake but will do the job either way")
        # if mirna_id:
        #     result = [self._mirna_retrieve(i, result) for i in mirna_id if self._mirna_retrieve(i, result) is not None]
        # if prec_id:
        #     for elem in prec_id:
        #         try:
        #             for mi in self._precursors_ID[elem].miRNAs:
        #                 if not self._exists(result, mi):
        #                     result.append(self._miRNAs_ID[mi])
        #         except:
        #             continue
        # if name:
        #     for mi in self._miRNAs_ID:
        #         if (self._miRNAs_ID[mi].mature_name == name) and not self._exists(result, self._miRNAs_ID[mi]):
        #             result.append(self._miRNAs_ID[mi])
        # if organism_name:
        #     for prec in self._precursors_ID:
        #         if self._precursors_ID[prec].organism == organism_name:
        #             mi_id = self._precursors_ID[prec].miRNAs
        #             for elem in mi_id:
        #                 if not self._exists(result, self._miRNAs_ID[elem]):
        #                     result.append(self._miRNAs_ID[elem])
        # if tax_level:
        #     for prec in self._precursors_ID:
        #         if tax_level in self._precursors_ID[prec].taxonomy:
        #             mi_id = self._precursors_ID[prec].miRNAs
        #             for elem in mi_id:
        #                 if not self._exists(result, self._miRNAs_ID[elem]):
        #                     result.append(self._miRNAs_ID[elem])
        # if chr:
        #     for mi in self._miRNAs_ID:
        #         if (chr in self._miRNAs_ID[mi].chromosome_mi) and not self._exists(result, self._miRNAs_ID[mi]):
        #             result.append(self._miRNAs_ID[mi])
        # if start and end:
        #     int_start = int(start)
        #     int_end = int(end)
        #     if not int_start < int_end:
        #         print(f"{Fore.RED}[Mir-Us]   Wrong coordinates; start cannot be lower than end")
        #         return None
        #     for mi in self._miRNAs_ID:
        #         for key in self._miRNAs_ID[mi].genome_coordinates_mi:
        #             for coord in self._miRNAs_ID[mi].genome_coordinates_mi[key]:
        #                 try:
        #                     f_start = int(coord[0])
        #                     f_stop = int(coord[1])
        #                 except:
        #                     continue
        #                 if (int_start <= f_start < int_end) and (int_start < f_stop <= int_end) and not self._exists(
        #                         result, self._miRNAs_ID[mi]):
        #                     result.append(self._miRNAs_ID[mi])
        # if strand:
        #     for mi in self._miRNAs_ID:
        #         if (strand in self._miRNAs_ID[mi].strand_mi) and not self._exists(result, self._miRNAs_ID[mi]):
        #             result.append(self._miRNAs_ID[mi])
        if not dict_result:
            return None
        if not len(dict_result) > 1:
            output = list(dict_result.values())[0]
            return output, len(output)
        print(f"{Fore.YELLOW}[Mir-Us]   Some of the given criteria were contradicting for the search system. Because of"
              f" that, the results are returned in a dictionary, where contradicting results are separated into"
              f" different search types. This search consist of (keys of generated dictionary): ")
        for key in dict_result.keys():
            print(f"{Fore.BLUE}'{key}': {len(dict_result[key])} results", sep="\n")
        #print([key for key in dict_result.keys()], sep="\n")
        return dict_result, sum([len(res) for res in list(dict_result.values())])

    @time_this
    def find_cluster(self, mirna_id=None, prec_id=None, search_type="up-downstream", range=None):
        """
        Returns all miRNAs present within given range from given miRNA in organism genome (dependent from given miRNA)
        :type mirna_id: str
        :param mirna_id: string representing miRNA's id
        :type prec_id: str
        :param prec_id: string representing precursor's id
        :type search_type: str
        :param search_type: string or int representing search type. The types are:
        - "up-downstream" - the search will be conducted below and above given miRNA position within given range
        ("up-downstream" is default)
        - "upstream" - the search will be conducted only above given miRNA position within given range
        - "downstream" - the search will be conducted only below given miRNA position within given range
        :type range: str
        :param range: string or int representing length of genome to be searched
        :return: list of miRNA objects
        :rtype: list
        """
        result = []

        def search_prec(self, start, org, range):
            int_start = 0
            int_end = 0
            count = 0
            if search_type == "up-downstream":
                int_start = start - range
                #print(int_start)
                int_end = start + range
                #print(int_end)
            if search_type == "upstream":
                int_start = start
                int_end = start + range
            if search_type == "downstream":
                int_start = start - range
                # print(int_start)
                int_end = start
            for prec in self._precursors_ID:
                if self._precursors_ID[prec].organism == org:
                    mi_id = self._precursors_ID[prec].miRNAs
                    for mi in mi_id:
                        for key in self._miRNAs_ID[mi].genome_coordinates_mi:
                            for coord in self._miRNAs_ID[mi].genome_coordinates_mi[key]:
                                try:
                                    f_start = int(coord[0])
                                    f_stop = int(coord[1])
                                except:
                                    continue
                                if (int_start <= f_start < int_end) and (
                                        int_start < f_stop <= int_end) and not self._exists(
                                        result, self._miRNAs_ID[mi]):
                                    count += 1
                                    print(f"{count}, {self._miRNAs_ID[mi].mature_ID}: {coord}")
                                    result.append(self._miRNAs_ID[mi])

        def search_prec2(self, start, org, range):
            int_start = 0
            int_end = 0
            count = 0
            if search_type == "up-downstream":
                int_start = start - range
                #print(int_start)
                int_end = start + range
                #print(int_end)
            if search_type == "upstream":
                int_start = start
                int_end = start + range
            if search_type == "downstream":
                int_start = start - range
                # print(int_start)
                int_end = start
            for prec in self._precursors_ID:
                if self._precursors_ID[prec].organism == org:
                    for coord in self._precursors_ID[prec].genome_coordinates:
                        try:
                            f_start = int(coord[0])
                            f_stop = int(coord[1])
                        except:
                            continue
                        if (int_start <= f_start < int_end) and (
                                int_start < f_stop <= int_end) and not self._exists(
                                result, self._precursors_ID[prec]):
                            count += 1
                            print(f"{count}, {self._precursors_ID[prec].precursor_ID}: {coord}")
                            result.append(self._precursors_ID[prec].precursor_ID)

        def search_mirna(self, start, org, range):
            int_start = 0
            int_end = 0
            count = 0
            if search_type == "up-downstream":
                int_start = start - range
                print(f"new_start: {int_start}")
                int_end = start + range
                print(f"new_end: {int_end}")
            if search_type == "upstream":
                int_start = start
                int_end = start + range
            if search_type == "downstream":
                int_start = start - range
                # print(int_start)
                int_end = start
            for mi in self._miRNAs_ID:
                if self._precursors_ID[self._miRNAs_ID[mi].precursor[0]].organism == org:
                    for key in self._miRNAs_ID[mi].genome_coordinates_mi:
                        for coord in self._miRNAs_ID[mi].genome_coordinates_mi[key]:
                            try:
                                #print(coord)
                                f_start = int(coord[0])
                                f_stop = int(coord[1])
                            except:
                                continue
                            if (int_start <= f_start < int_end) and (
                                    int_start < f_stop <= int_end) and not self._exists(
                                    result, self._miRNAs_ID[mi]):
                                count += 1
                                print(f"{count}, {self._miRNAs_ID[mi].mature_ID}: {coord}")
                                result.append(self._miRNAs_ID[mi])

        if mirna_id:
            org = self._precursors_ID[self._miRNAs_ID[mirna_id].precursor[0]].organism
            for key in self._miRNAs_ID[mirna_id].genome_coordinates_mi:
                for coord in self._miRNAs_ID[mirna_id].genome_coordinates_mi[key]:
                    start_position = int(coord[0])
                    print(f"org_start: {start_position}")
                    search_mirna(self, start_position, org, int(range))
        if prec_id:
            org = self._precursors_ID[prec_id].organism
            for coord in self._precursors_ID[prec_id].genome_coordinates:
                start_position = int(coord[0])
                #search_prec(self, start_position, org, int(range))
                search_prec2(self, start_position, org, int(range))
        if not result:
            return None

        return result, len(result)

    @time_this
    def get_tree(self):
        """
        Returns taxonomy tree, where each taxonomy level is a dictionary (nested dictionaries as access to another
        taxonomy level and list of organisms at particular taxonomy level if has any)
        :return: dictionary representing full taxonomy tree of organisms present in the base
        :rtype: dict
        """
        organism_codes = list(map(lambda x: getattr(x, "name"), self._organisms))  # lista wszystkich kodów organizmow
        tax_codes = list(map(lambda x: getattr(x, "tree"), self._organisms))
        tax_codes = [x.rstrip(";").split(";") for x in tax_codes]
        tax_dct = dict(zip(organism_codes, tax_codes))

        # structure def
        def tree():
            return dd(tree)

        # function to create tree (nested default dicts)
        def add_tree(t, keys):
            for key in keys:
                t = t[key]

        # function to access keys
        def get_tree(t, keys):
            for key in keys:
                t = t[key]
            return t

        # create tree with taxonomy
        org_tree = tree()
        for key in tax_dct:
            add_tree(org_tree, tax_dct[key])
        # org_tree = dd(list)
        # org_tree.default_factory = list()

        # add organisms
        for key in tax_dct:
            get_tree(org_tree, tax_dct[key])["!organism"] = [key]
        for key in tax_dct:
            if key not in get_tree(org_tree, tax_dct[key])["!organism"]:
                get_tree(org_tree, tax_dct[key])["!organism"].append(key)

        # sort in organism lists
        for key in tax_dct:
            get_tree(org_tree, tax_dct[key])["!organism"].sort()

        # debug
        print(json.dumps(org_tree, indent=4, sort_keys=True))
        return dict(org_tree), 1


class MiRLoad(MiRBase):
    """
    Parsing/loading class, which contains:
    - all functions to parse information directly from MiRBase using ftp
    - only newest version of MiRBase is supported
    """

    def __init__(self, version):
        super().__init__(version)
        #MiRBase.__init__(self, version)
        # self.base = MiRBase()

    def load_organisms(self, file_path):
        """
        Parses file with organisms (for 'current' database version - organisms.txt.gz)
        :param file_path: ftp path to files
        :type file_path: str
        :return: Creates namedtuples of organisms
        """
        with urllib.request.urlopen(file_path) as organisms_file_gz:
            with gzip.open(organisms_file_gz, mode='rt') as organisms_file:
                for line in organisms_file:
                    if not line.startswith("#"):
                        tmp = line.split('\t')
                        self._org_sh[tmp[0]] = tmp[2]
                        org = self._Organism(tmp[0], tmp[1], tmp[2], tmp[3], tmp[4].rstrip())
                        self._organisms.append(org)

    def load_miRNA(self, file_path):
        """
        Parses file with miRNA data (for 'current' database version - miRNA.dat.gz)
        :param file_path: ftp path to files
        :type file_path: str
        :return: Creates Precursor and miRNA objects
        """
        with urllib.request.urlopen(file_path) as miRNA_file_gz:
            with gzip.open(miRNA_file_gz, mode='rt') as miRNA_file:
                for line in miRNA_file:
                    if line.startswith("ID"):
                        name_p = line.split()[1]
                        while not line.startswith("AC"):
                            line = next(miRNA_file)
                        id_p = line.split()[1][:-1]
                        while not line.startswith("DE"):
                            line = next(miRNA_file)
                        org = name_p.split("-")[0]
                        full_name = self._org_sh[org]
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
                                products.append([start, end])
                            elif line.startswith("FT") and "/accession=" in line:
                                ac = line.split("=")[1].replace('\"', '').strip()
                                products[c].append(ac)
                            elif line.startswith("FT") and "/product=" in line:
                                name_mat = line.split("=")[1].replace('\"', '').strip()
                                products[c].append(name_mat)
                            elif line.startswith("FT") and "/evidence=" in line:
                                ev = line.split("=")[1].replace('\"', '').strip()
                                products[c].append(ev)
                                if ev == "experimental":
                                    line = next(miRNA_file)
                                    ex = ""
                                    if line.startswith("FT") and "/experiment=" in line:
                                        e = line.split("=")[1].replace('\"', '').strip()
                                        ex += e
                                        if line.strip()[-1] != "\"":
                                            line = next(miRNA_file)
                                            next_line = False
                                            while line.startswith(
                                                    "FT") and "/experiment=" not in line and "miRNA" not in line:
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
                        if id_p not in self._precursors_ID:
                            prec = miObject.Precursor(id_p, name_p, seq_full, full_name, ref, m_ids)
                            # prec = miRNA_draft.Precursor(id_p, name_p, seq_full, full_name, ref, m_ids).to_json()
                            # prec = miRNA_draft.Precursor(id_p, name_p, seq_full, full_name, ref, m_ids).dumper()
                            self._precursors_ID[id_p] = prec
                            self._precursors_name[name_p] = id_p
                        else:
                            for m in m_ids:
                                if m not in self._precursors_ID[id_p].miRNAs:
                                    self._precursors_ID[id_p].miRNAs.append(m)
                        for miRNA_vals in products:
                            if "5p" in miRNA_vals[3]:
                                e = "5p"
                            elif "3p" in miRNA_vals[3]:
                                e = "3p"
                            else:
                                e = "-"
                            if miRNA_vals[2] not in self._miRNAs_ID:
                                # mature = miRNA_draft.MiRNA([id_p],[miRNA_vals[2]],[miRNA_vals[3]],[miRNA_vals[0],miRNA_vals[1]],[miRNA_vals[4]],[miRNA_vals[5]],[e])
                                mature = miObject.MiRNA([id_p], [miRNA_vals[2]], [miRNA_vals[3]], full_name,
                                                        (miRNA_vals[0], miRNA_vals[1]), [miRNA_vals[4]],
                                                        [miRNA_vals[5]], [e], ref)
                                mature.get_mature_seq(seq_full, [miRNA_vals[0], miRNA_vals[1]])
                                self._miRNAs_ID[miRNA_vals[2]] = mature
                                # self._miRNAs_ID[miRNA_vals[2]] = mature.to_json()
                                # self._miRNAs_ID[miRNA_vals[2]] = mature.dumper()
                            else:
                                self._miRNAs_ID[miRNA_vals[2]].precursor.append(id_p)
                                self._miRNAs_ID[miRNA_vals[2]].mature_name.append(miRNA_vals[3])
                                self._miRNAs_ID[miRNA_vals[2]].mature_positions.append((miRNA_vals[0], miRNA_vals[1]))
                                self._miRNAs_ID[miRNA_vals[2]].evidence.append(miRNA_vals[4])
                                self._miRNAs_ID[miRNA_vals[2]].experiment.append(miRNA_vals[5])
                                self._miRNAs_ID[miRNA_vals[2]].end.append(e)
                                if len(ref) > 1:
                                    self._miRNAs_ID[miRNA_vals[2]].references.append(ref)
                                else:
                                    self._miRNAs_ID[miRNA_vals[2]].references.append(''.join(map(str, ref)))
                                # self._miRNAs_ID[miRNA_vals[2]].references = ref
                                self._miRNAs_ID[miRNA_vals[2]].get_mature_seq(seq_full, [miRNA_vals[0], miRNA_vals[1]])

    def load_hc(self, file_path):
        """
        Parses file with high-confidence data (for 'current' database version - hairpin_high_conf.fa.gz)
        :param file_path: ftp path to files
        :type file_path: str
        :return: Complements Precursor objects with information about high confidence of hairpin existence
        """
        with urllib.request.urlopen(file_path) as hc_file_gz:
            with gzip.open(hc_file_gz, mode='rt') as hc_file:
                for line1 in hc_file:
                    if line1.startswith('>'):
                        tmp_hc = line1.split()[1]
                        self._precursors_ID[tmp_hc].high_confidence = True
                        self._high_conf.append(tmp_hc)

    def load_structures(self, file_path):
        """
        Parses file with structure data (for 'current' database version - miRNA.str.gz)
        :param file_path: ftp path to files
        :type file_path: str
        :return: Complements Precursor objects with structures in dot-bracket format
        """
        # na podstawie pliku .str tworzy uzupelnia obiekty Precursor o strukture w notacji dot-bracket
        iupac = ['a', 'c', 'g', 't', 'u', 'r', 'y', 's', 'w', 'k', 'm', 'b', 'd', 'h', 'v', 'n']
        start = 0
        end = 7
        x = 0
        with urllib.request.urlopen(file_path) as struct_file_gz:
            struct_file = gzip.open(struct_file_gz, mode='rt')
            # print("THE FILE IS BEING PROCESSED...")
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
                    acgu_dt = ''.join(acgu_dt).lower().replace("-",
                                                               "")  # .replace('a', '.').replace('c', '.').replace('g','.').replace(
                    # 'u', '.').replace('n', '.').replace('t', '.').replace('r', '.').replace('y','.').replace('s', '.').replace(
                    # 'w', '.').replace('k', '.').replace('m', '.').replace('b', '.').replace('d', '.').replace('h', '.').replace(
                    # 'v', '.')
                    for nt in iupac:
                        if nt in acgu_dt:
                            acgu_dt = acgu_dt.replace(nt, '.')
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

                id_prec = self._precursors_name[name2]
                self._structures[id_prec] = all_dt_seq

            struct_file.close()

    def load_taxonomy(self):
        """
        Uses existing information to create dictionary of taxonomy
        :return: Complements Precursor objects with information about taxonony of organisms
        """
        organism_codes = list(map(lambda x: getattr(x, "name"), self._organisms))  # lista wszystkich kodów organizmow
        tax_codes = list(map(lambda x: getattr(x, "tree"), self._organisms))
        tax_dct = dict(zip(organism_codes, tax_codes))
        print(tax_dct)

        for prec in self._precursors_ID:
            p_id = self._precursors_ID[prec].precursor_ID
            p_org = self._precursors_ID[prec].organism
            self._organisms_of_prec[p_org].append(p_id)
            tax_name = tax_dct[p_org].split(';')[:-1]
            # print(tax_name)
            self._precursors_ID[prec].taxonomy = tax_name
            for tax_level in tax_name:
                if not (p_id in self._taxonomy_of_prec[tax_level]):
                    self._taxonomy_of_prec[tax_level].append(p_id)
            # self._precursors_ID[prec].structure = self._structures[prec]

    def load_genome(self, file_path):
        """
        Parses file with genomic data (for 'current' database version - /genomes/<organism>.gff3)
        :param file_path: ftp path to files
        :type file_path: str
        :return: Complements Precursor and miRNA objects with genomic information (chromosome, strand, coordinates)
        """
        for org in self._org_sh:
            try:
                with urllib.request.urlopen(file_path + org + '.gff3') as miRNA_file:
                    file2 = miRNA_file.readlines()
            except:
                continue
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
                    # sim_ID = split_info_miRNA[0].split('=')[1]
                    sim_Alias = split_info_miRNA[1].split('=')[1]
                    # sim_name = split_info_miRNA[2].split('=')[1]

                    if miRNA_type == 'miRNA_primary_transcript':
                        try:
                            self._precursors_ID[sim_Alias].chromosome.append(chr_miRNA)
                            self._precursors_ID[sim_Alias].strand.append(strand_seq)
                            self._precursors_ID[sim_Alias].genome_coordinates.append((start_seq, end_seq))
                        except:
                            continue
                    elif miRNA_type == 'miRNA':
                        try:
                            derivative = split_info_miRNA[3].split('=')[1].rstrip('\n')
                            self._miRNAs_ID[sim_Alias].chromosome_mi.append(chr_miRNA)
                            self._miRNAs_ID[sim_Alias].strand_mi.append(strand_seq)
                            # self._miRNAs_ID[sim_Alias].genome_coordinates_mi.append((start_seq, end_seq))
                            self._miRNAs_ID[sim_Alias].genome_coordinates_mi[derivative].append((start_seq, end_seq))
                        except:
                            continue
