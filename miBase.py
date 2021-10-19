"""
Main file of the package. All of the actions are handled with this code:
- parsing mirBase files and compiling them to .mir files
- reading data from .mir files
- accessing data from .mir file by package user
"""

import gzip
import json
import operator
import os
import urllib.request
from collections import defaultdict as dd
from collections import namedtuple as nt
from functools import reduce
from timeit import default_timer as timer

import dill
from Bio import Entrez
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from alive_progress import alive_bar
from colorama import init, Fore

import miObject
import utils

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


# CLASSES
class MiRBase:
    """Main operational class. This class defines which version of database will be used. Defines all available
    functions to retrieve data by user.

    """

    def __init__(self, version="CURRENT"):
        """Initializes MiRBase object from which data can be accessed using provided functions. Database version might
        be specified.

        !!! info "There are multiple versions of miRBase database compatible with Mir-Us. [Details :octicons-link-16:](versions.md){: target="_blank" .md-button .md-button--primary }"

        Args:
            version (str): ID of database version. "CURRENT" is the most recent version (22.1).
        """
        self._ftp_path = "ftp://mirbase.org/pub/mirbase/"  # main path to all files in mirbase ftp
        self._miRBase_version = version  # version on which current instance of tool will be working
        self._miRNAs_ID = {}  # dict of miRNA ID : miRNA object
        self._precursors_ID = {}  # dict of pre-miRNA_ID : Precursor object
        self._org_sh = {}  # dict of 3-letter code : organism full name
        self._organisms = []  # list of namedtuples Organism = nt('Organism', 'organism division name tree taxid')
        self._taxonomy_of_prec = dd(list)  # dict of every tax level with list of precursor and/or miRNAs IDs
        self._organisms_of_prec = dd(list)  # dict of every organism with list of precursors and/or miRNAs IDs
        self._high_conf = []  # list of high confidence pre-miRNAs and/or miRNAs (IDs)
        self._structures = {}  # dict of pre-miRNA ID and their structures from miRNA.str file
        self._precursors_name = {}  # dict of pre-miRNA name : pre-miRNA ID
        self._matures_name = {}  # dict of miRNA name : miRNA ID

        self._versions = utils._cache_versions()

        self._loader = MiRLoad  # reference to loader class - very important

        self._Organism = nt('Organism', 'organism division name tree taxid')

        utils._show_banner()

        # check files integrity
        try:
            os.makedirs(f"{os.getcwd()}/data", exist_ok=True)
            for elem in self._versions.keys():
                path = os.path.join("data", elem)
                os.makedirs(path, exist_ok=True)
            self._load_all_data()
        except:
            print(f"{Fore.RED}[Mir-Us]   Missing data files (.mir); "
                  f"{Fore.YELLOW}performing compilation of data files...\n"
                  f"{Fore.YELLOW}Please, be patient, compiling will take several minutes.")
            try:
                self._compile_indexes()
                init()
                print(f"{Fore.YELLOW}[Mir-Us]   Data files compiled successfully!")
                self._load_all_data()
            except Exception as e:
                init()
                print(
                    utils._fatal_error_handle(
                        e) + "Log file with full error description is located in /logs directory.")
                exit()

    def _compile_indexes(self):
        """Produces .mir files which contain indexed data.
        """
        path = self._ftp_path + self._miRBase_version

        # INITIALIZE PROGRESSBAR-------------------------------------
        with alive_bar(10, bar='blocks') as bar:
            # COMPILE ORGANISMS--------------------------------------
            self._loader.load_organisms(self, file_path=path + self._versions[self._miRBase_version]["org_file"])
            with open(f'data/{self._miRBase_version}/organisms.mir', 'wb') as fh_org_dump:
                dill.dump(self._organisms, fh_org_dump)
            fh_org_dump.close()
            bar()
            with open(f'data/{self._miRBase_version}/org_short.mir', 'wb') as fh_orgsh_dump:
                dill.dump(self._org_sh, fh_orgsh_dump)
            fh_orgsh_dump.close()
            bar()
            # -------------------------------------------------------

            # COMPILE MIRNA------------------------------------------
            self._loader.load_miRNA(self, file_path=path + self._versions[self._miRBase_version]["mirna_dat"])
            self._loader.load_genome(self, file_path=path + self._versions[self._miRBase_version]["genomes"])
            with open(f'data/{self._miRBase_version}/precursors_ID.mir', 'wb') as fh_precid_dump:
                dill.dump(self._precursors_ID, fh_precid_dump)
            fh_precid_dump.close()
            bar()
            with open(f'data/{self._miRBase_version}/precursors_name.mir', 'wb') as fh_precname_dump:
                dill.dump(self._precursors_name, fh_precname_dump)
            fh_precname_dump.close()
            bar()
            with open(f'data/{self._miRBase_version}/miRNAs_ID.mir', 'wb') as fh_mirnasid_dump:
                dill.dump(self._miRNAs_ID, fh_mirnasid_dump)
            fh_mirnasid_dump.close()
            bar()
            with open(f'data/{self._miRBase_version}/matures_name.mir', 'wb') as fh_maturename_dump:
                dill.dump(self._matures_name, fh_maturename_dump)
            fh_maturename_dump.close()
            bar()
            # -------------------------------------------------------

            # COMPILE HIGH-CONF--------------------------------------
            try:
                self._loader.load_hc(self, file_path=path + self._versions[self._miRBase_version]["high_conf"])
                with open(f'data/{self._miRBase_version}/high_conf.mir', 'wb') as fh_high_dump:
                    dill.dump(self._high_conf, fh_high_dump)
                fh_high_dump.close()
            except:
                pass
            bar()
            # -------------------------------------------------------

            # COMPILE STRUCTURES-------------------------------------
            self._loader.load_structures(self, file_path=path + self._versions[self._miRBase_version]["mirna_str"])
            with open(f'data/{self._miRBase_version}/structures.mir', 'wb') as fh_mirstruc_dump:
                dill.dump(self._structures, fh_mirstruc_dump)
            fh_mirstruc_dump.close()
            bar()
            # -------------------------------------------------------

            # COMPILE TAXONOMY---------------------------------------
            self._loader.load_taxonomy(self)
            with open(f'data/{self._miRBase_version}/taxonomy_prec.mir', 'wb') as fh_taxprec_dump:
                dill.dump(self._taxonomy_of_prec, fh_taxprec_dump)
            fh_taxprec_dump.close()
            bar()
            with open(f'data/{self._miRBase_version}/taxonomy_org.mir', 'wb') as fh_taxorg_dump:
                dill.dump(self._organisms_of_prec, fh_taxorg_dump)
            fh_taxorg_dump.close()
            bar()
            # -------------------------------------------------------

    def _load_all_data(self):
        """Loads all data from .mir files and assigns to data structures
        """
        # LOAD ORGANISMS-----------------------------------------
        with open(f'data/{self._miRBase_version}/organisms.mir', 'rb') as fh_org_load:
            self._organisms = dill.load(fh_org_load)
        fh_org_load.close()
        # -------------------------------------------------------

        # LOAD MIRNA---------------------------------------------
        with open(f'data/{self._miRBase_version}/precursors_ID.mir', 'rb') as fh_precid_load:
            self._precursors_ID = dill.load(fh_precid_load)
        fh_precid_load.close()

        with open(f'data/{self._miRBase_version}/precursors_name.mir', 'rb') as fh_precname_load:
            self._precursors_name = dill.load(fh_precname_load)
        fh_precid_load.close()

        with open(f'data/{self._miRBase_version}/miRNAs_ID.mir', 'rb') as fh_mirnasid_load:
            self._miRNAs_ID = dill.load(fh_mirnasid_load)
        fh_mirnasid_load.close()

        with open(f'data/{self._miRBase_version}/matures_name.mir', 'rb') as fh_maturename_load:
            self._matures_name = dill.load(fh_maturename_load)
        fh_maturename_load.close()

        with open(f'data/{self._miRBase_version}/org_short.mir', 'rb') as fh_orgsh_load:
            self._org_sh = dill.load(fh_orgsh_load)
        fh_orgsh_load.close()
        # -------------------------------------------------------

        # LOAD HIGH-CONF-----------------------------------------
        if self._versions[self._miRBase_version]["high_conf"] is not None:
            with open(f'data/{self._miRBase_version}/high_conf.mir', 'rb') as fh_high_load:
                self._high_conf = dill.load(fh_high_load)
            fh_high_load.close()
        # -------------------------------------------------------

        # LOAD STRUCTURES----------------------------------------
        with open(f'data/{self._miRBase_version}/structures.mir', 'rb') as fh_mirstruc_load:
            self._structures = dill.load(fh_mirstruc_load)
        fh_mirstruc_load.close()
        # -------------------------------------------------------

        # LOAD TAXONOMY------------------------------------------
        with open(f'data/{self._miRBase_version}/taxonomy_prec.mir', 'rb') as fh_taxprec_load:
            self._taxonomy_of_prec = dill.load(fh_taxprec_load)
        fh_taxprec_load.close()
        with open(f'data/{self._miRBase_version}/taxonomy_org.mir', 'rb') as fh_taxorg_load:
            self._organisms_of_prec = dill.load(fh_taxorg_load)
        fh_taxorg_load.close()
        # -------------------------------------------------------

        # PERFORM MERGE------------------------------------------
        self._merge_data()
        # -------------------------------------------------------

    def _merge_data(self):
        """Merges data - some data structures lacks certain parts and this function completes them.
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
            if not utils._exists(current_result, self._precursors_ID[passed_id]):
                return self._precursors_ID[passed_id]
        except KeyError as e:
            pass

    def _mirna_retrieve(self, passed_id, current_result):
        try:
            if not utils._exists(current_result, self._miRNAs_ID[passed_id]):
                return self._miRNAs_ID[passed_id]
        except KeyError as e:
            pass

    def _make_tax_dict(self):
        organism_codes = list(map(lambda x: getattr(x, "name"), self._organisms))
        tax_codes = list(map(lambda x: getattr(x, "tree"), self._organisms))
        tax_codes = [x.rstrip(";").split(";") for x in tax_codes]
        tax_dct = dict(zip(organism_codes, tax_codes))
        return tax_dct

    def get_organisms_list(self):
        """Returns list of all organisms.

        Returns:
            list[namedtuple]: Organism namedtuple which is organized as follows: `Organism(organism='', division='', name='', tree='', taxid='')`

        ???+ info "Organism namedtuple structure explanation"

            | Key          | Explanation                 |
            | ------------ | --------------------------- |
            | `organism`   | Organism name abbreviation. |
            | `division`   | Organism name abbreviation. |
            | `name`       | Full organism name in latin. |
            | `tree`       | Full taxonomy path of the organism delimited by `;` |
            | `taxid`      | NCBI taxonomy ID of the organism |
        """
        return self._organisms

    @utils.time_this
    def get_tax_level(self, organism=None, verbose=False):
        """Returns taxonomy level assigned to organism.

        Args:
            organism (list[str]): list of organisms names
            verbose (bool): A flag, which allows or disallows showing search details (number of returned elements, time of execution, errors)

        Returns:
            Optional[dict]: Dictionary of organisms and its assigned taxonomy, which is a list of tax levels
            (key: organism, value: taxonomy) or `None` if no results are found.
        """
        if isinstance(organism, str):
            organism = [organism]
        result = {}
        tax_dct = self._make_tax_dict()
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

    @utils.time_this
    def get_organism(self, tax=None, verbose=False):
        """Returns organisms which are assigned to a given taxonomy level.

        Args:
            tax (str): Full name of taxonomy level
            verbose (bool): A flag, which allows or disallows showing search details (number of returned elements, time of execution, errors)

        Returns:
            Optional[list[str]]: Full names of organisms representing given taxonomy level or `None` if no results are
            found.
        """
        result = []
        tax_dct = self._make_tax_dict()
        # print(tax_dct)
        for key, value in tax_dct.items():
            if tax in value:
                result.append(key)
            else:
                continue
        if not result:
            return None
        return result, len(result)

    @utils.time_this
    def get_organisms_short(self, organism=None, verbose=False):
        """Returns dictionary of all organisms names and their abbreviations or a single organism abbreviation.

        Args:
            organism (str): Full organism name

        Returns:
            Optional[str, dict]: `str` which is an abbreviation of the given organism name or `dict` of all
            abbreviations and its full organism names with the following dictionary structure: (key: abbreviation,
            value: full organism name)
        """
        result = []
        if organism is not None:
            result.append(list(self._org_sh.keys())[list(self._org_sh.values()).index(organism)])
        else:
            result = self._org_sh
        if not result:
            return None
        return result[0] if organism else result, len(result)

    @utils.time_this
    def get_taxid(self, organism=None, verbose=False):
        """Returns NCBI taxonomy ID (taxid) assigned to organism

        Args:
            organism (list[str]): Full organism name
            verbose (bool): A flag, which allows or disallows showing search details (number of returned elements, time of execution, errors)

        Returns:
            Optional[dict]: Dictionary of organisms with assigned taxid (key: full organism name, value: taxid) or
            `None` if no results are found.
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

    @utils.time_this
    def get_precursor(self, prec_id: list = None, name="", organism_name="", tax_level="", chr="", start="",
                      end="", strand='', mirna_id=None, verbose=False):
        """Returns precursor objects according to a given search criteria

        Args:
            prec_id (list[str]): Precursor IDs
            name (str): Full precursor name
            organism_name (str): Full organism name in which precursor is present
            tax_level (str): Taxonomy level affiliated with precursor
            chr (str): Chromosome name in which precursor is present
            start (str): Genomic location in which precursor sequence starts
            end (str): Genomic location in which precursor sequence ends
            strand (str): Strand name in which precursor is present
            mirna_id (list[str]): Affiliated precursor miRNA IDs
            verbose (bool): A flag, which allows or disallows showing search details (number of returned elements, time of execution, errors)

        Returns:
            Optional[list[Precursor], dict]: `list` of Precursor objects matching given criteria **only if single
            search is conducted.** `dict` **only if multiple searches are conducted.** For each search type
            ([**additional information**](#additional-information)) a `list` of Precursor objects that matches criteria
            is returned. Eventually, `None` if no results are found.

        ## Additional information:
        ???+ info "Precursor search types"
            There are a few search 'types', depending on search criteria. That means, if the user choose a set of
            criteria that will be interpreted as contradicting, a separate search will be conducted for each 'type'.
            This is important, because multiple searches in single function call will return a dictionary with
            categorized search results.

            **Search types:**

            | Search type      | Set of criteria                                         |
            | ---------------- | ------------------------------------------------------- |
            | `mirna_id-search`| Only `mirna_id`                                         |
            | `prec_id-search` | Only `prec_id`                                          |
            | `name-search`    | Only `name`                                             |
            | `organism-search`| Only `organism_name`                                    |
            | `tax-search`     | Only `tax_level`                                        |
            | `chr-search`     | Only `chr`                                              |
            | `strand-search`  | Only `strand`                                           |
            | `genomic-search` | <ul><li>`organism`, `chr`</li><ul><li>`organism`, `chr`, `start`</li><li>`organism`, `chr`, `end`</li><li>`organism`, `chr`, `start`, `end`</li></ul><li>`organism`, `strand`</li><ul><li>`organism`, `strand`, `start`</li><li>`organism`, `strand`, `end`</li><li>`organism`, `strand`, `start`, `end`</li></ul><li>`organism`, `chr`, `strand`</li><ul><li>`organism`, `chr`, `strand`, `start`</li><li>`organism`, `chr`, `strand`, `end`</li><li>`organism`, `chr`, `strand`, `start`, `end`</li></ul></ul>|

            !!! warning "Search type is also a key in returned dictionary which allows to access found data."

            **Possible function returns:**

            | Condition                                 | Returned structure                   |
            | ----------------------------------------- | ------------------------------------ |
            | Single search (only one type of search)   | `list` of found objects or `None`    |
            | Multiple searches (contradicting types)   | `dict` of 'type: [found objects]' or `None`. If certain type of search was unsuccessful, its value will be `None` |
        """

        if isinstance(prec_id, str) and prec_id:
            prec_id = [prec_id]
        if isinstance(mirna_id, str) and mirna_id:
            mirna_id = [mirna_id]
        dict_result = {}
        first = True

        #if not end and start or start and end or organism_name and chr or organism_name and strand:
        if organism_name and chr or organism_name and strand:
            result = []
            if organism_name:
                for prec in self._precursors_ID:
                    if (self._precursors_ID[prec].organism == organism_name) and not \
                            utils._exists(result, self._precursors_ID[prec]):
                        result.append(self._precursors_ID[prec])
                first = False
            if chr:
                temp_result = []
                if first:
                    for prec in self._precursors_ID:
                        if (chr in self._precursors_ID[prec].chromosome) and not \
                                utils._exists(result, self._precursors_ID[prec]):
                            result.append(self._precursors_ID[prec])
                    first = False
                elif not first:
                    for res in result:
                        if chr in res.chromosome and not utils._exists(temp_result, res):
                            temp_result.append(res)
                result = temp_result
            if strand:
                if not utils._check_strand(strand):
                    print(f"{Fore.RED}[Mir-Us]   Incorrect 'strand' value; 'strand' can only be '+' or '-'")
                    return None
                temp_result = []
                if first:
                    for prec in self._precursors_ID:
                        if (strand in self._precursors_ID[prec].strand) and not utils._exists(result,
                                                                                              self._precursors_ID[
                                                                                                  prec]):
                            result.append(self._precursors_ID[prec])
                    first = False
                elif not first:
                    for res in result:
                        if strand in res.strand and not utils._exists(temp_result, res):
                            temp_result.append(res)
                result = temp_result
            if start and end:
                temp_result = []
                if first:
                    int_start = int(start)
                    int_end = int(end)
                    if int_start < 0 or int_end < 0:
                        print(f"{Fore.RED}[Mir-Us]   Incorrect 'start' or 'end' value; 'start' or 'end' cannot be less than zero.")
                        return None
                    if not int_start < int_end:
                        print(f"{Fore.RED}[Mir-Us]   Wrong coordinates; start cannot be lower than end")
                        return None
                    for prec in self._precursors_ID:
                        for coord in self._precursors_ID[prec].genome_coordinates:
                            try:
                                f_start = int(coord[0])
                                f_stop = int(coord[1])
                            except:
                                continue
                            if (int_start <= f_start < int_end) and (int_start < f_stop <= int_end) \
                                    and not utils._exists(result, self._precursors_ID[prec]):
                                result.append(self._precursors_ID[prec])
                    first = False
                elif not first:
                    int_start = int(start)
                    int_end = int(end)
                    if int_start < 0 or int_end < 0:
                        print(f"{Fore.RED}[Mir-Us]   Incorrect 'start' or 'end' value; 'start' or 'end' cannot be less than zero.")
                        return None
                    if not int_start < int_end:
                        print(f"{Fore.RED}[Mir-Us]   Wrong coordinates; start cannot be lower than end")
                        return None
                    for res in result:
                        for coord in res.genome_coordinates:
                            try:
                                f_start = int(coord[0])
                                f_stop = int(coord[1])
                            except:
                                continue
                            if (int_start <= f_start < int_end) and (int_start < f_stop <= int_end) \
                                    and not utils._exists(temp_result, res):
                                temp_result.append(res)
                    result = temp_result
            if not end and start:
                temp_result = []
                if first:
                    int_start = int(start)
                    if int_start < 0:
                        print(f"{Fore.RED}[Mir-Us]   Incorrect 'start' or 'end' value; 'start' or 'end' cannot be less than zero.")
                        return None
                    for prec in self._precursors_ID:
                        for coord in self._precursors_ID[prec].genome_coordinates:
                            try:
                                f_start = int(coord[0])
                            except:
                                continue
                            if (int_start <= f_start) and not utils._exists(result, self._precursors_ID[prec]):
                                result.append(self._precursors_ID[prec])
                    first = False
                elif not first:
                    int_start = int(start)
                    if int_start < 0:
                        print(f"{Fore.RED}[Mir-Us]   Incorrect 'start' or 'end' value; 'start' or 'end' cannot be less than zero.")
                        return None
                    for res in result:
                        for coord in res.genome_coordinates:
                            try:
                                f_start = int(coord[0])
                            except:
                                continue
                            if int_start <= f_start and not utils._exists(temp_result, res):
                                temp_result.append(res)
                    result = temp_result
            if not start and end:
                temp_result = []
                if first:
                    int_end = int(end)
                    if int_end < 0:
                        print(f"{Fore.RED}[Mir-Us]   Incorrect 'start' or 'end' value; 'start' or 'end' cannot be less than zero.")
                        return None
                    for prec in self._precursors_ID:
                        for coord in self._precursors_ID[prec].genome_coordinates:
                            try:
                                f_stop = int(coord[1])
                            except:
                                continue
                            if (int_end >= f_stop) and not utils._exists(result, self._precursors_ID[prec]):
                                result.append(self._precursors_ID[prec])
                    first = False
                elif not first:
                    int_end = int(end)
                    if int_end < 0:
                        print(f"{Fore.RED}[Mir-Us]   Incorrect 'start' or 'end' value; 'start' or 'end' cannot be less than zero.")
                        return None
                    for res in result:
                        for coord in res.genome_coordinates:
                            try:
                                f_stop = int(coord[1])
                            except:
                                continue
                            if int_end >= f_stop and not utils._exists(temp_result, res):
                                temp_result.append(res)
                    result = temp_result
            dict_result["genomic-search"] = result
        if prec_id:
            result = []
            result = [self._precursor_retrieve(i, result) for i in prec_id if
                      self._precursor_retrieve(i, result) is not None]
            dict_result["prec_id-serch"] = result
        if mirna_id:
            result = []
            for elem in mirna_id:
                try:
                    for prec in self._miRNAs_ID[elem].precursors:
                        if not utils._exists(result, prec):
                            result.append(self._precursors_ID[prec])
                except:
                    continue
            dict_result["mirna_id-search"] = result
        if name:
            result = []
            try:
                if not utils._exists(result, self._precursors_ID[self._precursors_name[name]]):
                    result.append(self._precursors_ID[self._precursors_name[name]])
            except:
                pass
            dict_result["name-search"] = result
        if not start and not end and not chr and not strand and organism_name:
            result = []
            for prec in self._precursors_ID:
                if (self._precursors_ID[prec].organism == organism_name) and not utils._exists(result,
                                                                                               self._precursors_ID[
                                                                                                   prec]):
                    result.append(self._precursors_ID[prec])
            dict_result["organism-search"] = result
        if tax_level:
            result = []
            for prec in self._precursors_ID:
                if (tax_level in self._precursors_ID[prec].taxonomy) and not utils._exists(result,
                                                                                           self._precursors_ID[prec]):
                    result.append(self._precursors_ID[prec])
            dict_result["tax-search"] = result
        if not start and not end and not organism_name and not strand and chr:
            result = []
            for prec in self._precursors_ID:
                if (chr in self._precursors_ID[prec].chromosome) and not utils._exists(result,
                                                                                       self._precursors_ID[prec]):
                    result.append(self._precursors_ID[prec])
            dict_result["chr-search"] = result
        if not start and not end and not organism_name and not chr and strand:
            if not utils._check_strand(strand):
                print(f"{Fore.RED}[Mir-Us]   Incorrect 'strand' value; 'strand' can only be '+' or '-'")
                return None
            result = []
            for prec in self._precursors_ID:
                if (strand in self._precursors_ID[prec].strand) and not utils._exists(result,
                                                                                      self._precursors_ID[prec]):
                    result.append(self._precursors_ID[prec])
            dict_result["strand-search"] = result
        if len(dict_result) >= 1 and not bool([res for res in dict_result.values() if res != []]):
            return None
        if not len(dict_result) > 1:
            # if list(dict_result.values())[0]:
            #if bool(dict_result):
            if dict_result.values():
                output = list(dict_result.values())[0]
                return output, len(output)
            else:
                return None
        # if len(dict_result) > 1 and not bool([res for res in dict_result.values() if res != []]):
        #     return None
        print(f"{Fore.YELLOW}[Mir-Us]   Some of the given criteria were contradicting for the search system. Because of"
              f" that, the results are returned in a dictionary, where contradicting results are separated into"
              f" different search types. This search consist of (keys of generated dictionary): ")
        for key in dict_result.keys():
            print(f"{Fore.BLUE}'{key}': {len(dict_result[key])} results", sep="\n")
        # print([key for key in dict_result.keys()], sep="\n")
        return dict_result, sum([len(res) for res in list(dict_result.values())])

    @utils.time_this
    def get_references(self, mirna_id=None, mirna_name=None, prec_id=None, prec_name=None, link=False, verbose=False):
        """Returns list of references

        Args:
            mirna_id (list[str]): IDs of miRNA
            mirna_name (list[str]): Full miRNA names
            prec_id (list[str]): IDs of precursors
            prec_name (list[str]): Full precursor names
            link (bool): A flag, which forces function to return PubMed links instead of reference numbers
            verbose (bool): A flag, which allows or disallows showing search details (number of returned elements, time of execution, errors)

        Returns:
            Optional[dict]: Dictionary of precursors or miRNAs ids or names with assigned lists of reference accession numbers or
             Pubmed links (depending on `link` flag) (key: precursor or miRNA name or id, value: list of accession numbers or
             links) or `None` if no results are found.
        """

        if isinstance(mirna_id, str):
            mirna_id = [mirna_id]
        if isinstance(mirna_name, str):
            mirna_name = [mirna_name]
        if isinstance(prec_id, str):
            prec_id = [prec_id]
        if isinstance(prec_name, str):
            prec_name = [prec_name]
        #result = []
        result = dd(list)
        if mirna_id:
            for m_id in mirna_id:
                try:
                    for ref in self._miRNAs_ID[m_id].references:
                        if ref not in result and link is not True:
                            #result.append(ref)
                            result[m_id].append(ref)
                        elif isinstance(link, bool) and link is True:
                            #result.append(f"https://pubmed.ncbi.nlm.nih.gov/{ref}/")
                            result[m_id].append(f"https://pubmed.ncbi.nlm.nih.gov/{ref}/")
                except:
                    continue
        if mirna_name:
            for mi_name in mirna_name:
                try:
                    for ref in self._miRNAs_ID[self._matures_name[mi_name]].references:
                        if ref not in result and link is not True:
                            #result.append(ref)
                            result[mi_name].append(ref)
                        elif isinstance(link, bool) and link is True:
                            #result.append(f"https://pubmed.ncbi.nlm.nih.gov/{ref}/")
                            result[mi_name].append(f"https://pubmed.ncbi.nlm.nih.gov/{ref}/")
                except:
                    continue
        if prec_name:
            for p_name in prec_name:
                try:
                    for ref in self._precursors_ID[self._precursors_name[p_name]].references:
                        if ref not in result and link is not True:
                            #result.append(ref)
                            result[p_name].append(ref)
                        elif isinstance(link, bool) and link is True:
                            #result.append(f"https://pubmed.ncbi.nlm.nih.gov/{ref}/")
                            result[p_name].append(f"https://pubmed.ncbi.nlm.nih.gov/{ref}/")
                except:
                    continue
            # for mi_name in mirna_name:
            #     try:
            #         for mi in self._miRNAs_ID:
            #             if mi_name in self._miRNAs_ID[mi].name:
            #                 for ref in self._miRNAs_ID[mi].references:
            #                     if ref not in result and link is not True:
            #                         result.append(ref)
            #                     elif isinstance(link, bool) and link is True:
            #                         result.append(f"https://pubmed.ncbi.nlm.nih.gov/{ref}/")
            #     except:
            #         continue
        if prec_id:
            for p_id in prec_id:
                try:
                    for ref in self._precursors_ID[p_id].references:
                        if ref not in result and link is not True:
                            #result.append(ref)
                            result[p_id].append(ref)
                        elif isinstance(link, bool) and link is True:
                            #result.append(f"https://pubmed.ncbi.nlm.nih.gov/{ref}/")
                            result[p_id].append(f"https://pubmed.ncbi.nlm.nih.gov/{ref}/")
                except:
                    continue
        if not result:
            return None
        return dict(result), len(result)

    @utils.time_this
    def get_structure(self, id=None, name=None, verbose=False):
        """Returns dictionary of precursors IDs with assigned structures in dot-bracket format

        Args:
            id (list[str]): Precursor IDs
            name (list[str]): Full precursor names
            verbose (bool): A flag, which allows or disallows showing search details (number of returned elements, time of execution, errors)

        Returns:
            Optional[dict]: Dictionary of precursors ids or names with assigned structure (key: precursor name or id, value: structure)
            or `None` if no results are found.
        """

        result = {}
        id_result = {}
        name_result = {}
        if isinstance(id, str):
            id = [id]
        if isinstance(name, str):
            name = [name]
        if id:
            try:
                # for i in id:
                #     try:
                #         id_result[self._precursors_ID[i].ID] = self._precursors_ID[i].structure
                #     except:
                #         continue
                id_result = {self._precursors_ID[i].ID: self._precursors_ID[i].structure for i in id}
            except:
                pass
        if name:
            try:
                #for i in name:
                    #try:
                name_result = {self._precursors_ID[prec].name: self._precursors_ID[prec].structure for
                          prec in self._precursors_ID if self._precursors_ID[prec].name in name}
                #name_result = dict_res
                    #except:
                    #    continue
            except:
                pass
        result = {**id_result, **name_result}
        if not result:
            return None
        return result, len(result)

    @utils.time_this
    def get_mirna(self, mirna_id: list = None, name="", organism_name="", tax_level="", chr="", start="",
                  end="", strand='', prec_id: list = None, verbose=False):
        """Returns miRNA objects according to given search criteria.

        Args:
            mirna_id (list[str]): miRNA IDs
            name (str): Full miRNA name
            organism_name (str): Full organism name affiliated with miRNA
            tax_level (str): Taxonomy level affiliated with miRNA
            chr (str):  Chromosome name in which miRNA is present
            start (str): Genomic location in which miRNA sequence starts
            end (str): Genomic location in which miRNA sequence ends
            strand (str): Strand name in which miRNA is present
            prec_id (list[str]): Affiliated precursor IDs
            verbose (bool): A flag, which allows or disallows showing search details (number of returned elements, time of execution, errors)

        Returns:
            Optional[list[MiRNA], dict]: `list` of miRNA objects matching criteria **only if single search is conducted.** `dict` **only if
            multiple searches are conducted.** For each search type ([**additional information**](#additional-information)) a `list` of miRNA objects that matches criteria is
             returned. Eventually, `None` if no results are found.

        ## Additional information:
        ???+ info "miRNA search types"
            There are a few search 'types', depending on search criteria. That means, if the user choose a set of
            criteria that will be interpreted as contradicting, a separate search will be conducted for each 'type'.
            This is important, because multiple searches in single function call will return a dictionary with
            categorized search results.

            **Search types:**

            | Search type      | Set of criteria                                         |
            | ---------------- | ------------------------------------------------------- |
            | `mirna_id-search`| Only `mirna_id`                                         |
            | `prec_id-search` | Only `prec_id`                                          |
            | `name-search`    | Only `name`                                             |
            | `organism-search`| Only `organism_name`                                    |
            | `tax-search`     | Only `tax_level`                                        |
            | `chr-search`     | Only `chr`                                              |
            | `strand-search`  | Only `strand`                                           |
            | `genomic-search` | <ul><li>`organism`, `chr`</li><ul><li>`organism`, `chr`, `start`</li><li>`organism`, `chr`, `end`</li><li>`organism`, `chr`, `start`, `end`</li></ul><li>`organism`, `strand`</li><ul><li>`organism`, `strand`, `start`</li><li>`organism`, `strand`, `end`</li><li>`organism`, `strand`, `start`, `end`</li></ul><li>`organism`, `chr`, `strand`</li><ul><li>`organism`, `chr`, `strand`, `start`</li><li>`organism`, `chr`, `strand`, `end`</li><li>`organism`, `chr`, `strand`, `start`, `end`</li></ul></ul>|

            !!! warning "Search type is also a key in returned dictionary which allows to access found data."

            **Possible function returns:**

            | Condition                                 | Returned structure                   |
            | ----------------------------------------- | ------------------------------------ |
            | Single search (only one type of search)   | `list` of found objects or `None`    |
            | Multiple searches (contradicting types)   | `dict` of 'type: [found objects]' or `None`. If certain type of search was unsuccessful, its value will be `None` |

        """

        if isinstance(mirna_id, str) and mirna_id:
            mirna_id = [mirna_id]
        if isinstance(prec_id, str) and prec_id:
            prec_id = [prec_id]
        # print(mirna_id)
        dict_result = {}
        first = True

        #if not end and start or start and end or organism_name and chr or organism_name and strand:
        if organism_name and chr or organism_name and strand:
            result = []
            if organism_name:
                for mi in self._miRNAs_ID:
                    if (self._miRNAs_ID[mi].organism == organism_name) and not utils._exists(result,
                                                                                             self._miRNAs_ID[mi]):
                        result.append(self._miRNAs_ID[mi])
                first = False
            if chr:
                temp_result = []
                if first:
                    for mi in self._miRNAs_ID:
                        if (chr in self._miRNAs_ID[mi].chromosome) and not utils._exists(result,
                                                                                            self._miRNAs_ID[mi]):
                            result.append(self._miRNAs_ID[mi])
                    first = False
                elif not first:
                    # count = 0
                    # print(len(result))
                    for res in result:
                        # count += 1
                        # print(f"Current {(count)}: {res.chromosome}")
                        if chr in res.chromosome and not utils._exists(temp_result, res):
                            # print(f"Kept:{res.chromosome}")
                            temp_result.append(res)
                result = temp_result
                # print("Chr from genomic context")
            if strand:
                if not utils._check_strand(strand):
                    print(f"{Fore.RED}[Mir-Us]   Incorrect 'strand' value; 'strand' can only be '+' or '-'")
                    return None
                temp_result = []
                if first:
                    for mi in self._miRNAs_ID:
                        if (strand in self._miRNAs_ID[mi].strand) and not utils._exists(result, self._miRNAs_ID[mi]):
                            result.append(self._miRNAs_ID[mi])
                    first = False
                elif not first:
                    for res in result:
                        if strand in res.strand and not utils._exists(temp_result, res):
                            temp_result.append(res)
                result = temp_result
                # print("Strand from genomic context")
            if start and end:
                temp_result = []
                if first:
                    int_start = int(start)
                    int_end = int(end)
                    if int_start < 0 or int_end < 0:
                        print(f"{Fore.RED}[Mir-Us]   Incorrect 'start' or 'end' value; 'start' or 'end' cannot be less than zero.")
                        return None
                    if not int_start < int_end:
                        print(f"{Fore.RED}[Mir-Us]   Wrong coordinates; start cannot be lower than end")
                        return None
                    for mi in self._miRNAs_ID:
                        for key in self._miRNAs_ID[mi].genome_coordinates:
                            for coord in self._miRNAs_ID[mi].genome_coordinates[key]:
                                try:
                                    f_start = int(coord[0])
                                    f_stop = int(coord[1])
                                except:
                                    continue
                                if (int_start <= f_start < int_end) and (int_start < f_stop <= int_end) \
                                        and not utils._exists(result, self._miRNAs_ID[mi]):
                                    result.append(self._miRNAs_ID[mi])
                    first = False
                elif not first:
                    int_start = int(start)
                    int_end = int(end)
                    if int_start < 0 or int_end < 0:
                        print(f"{Fore.RED}[Mir-Us]   Incorrect 'start' or 'end' value; 'start' or 'end' cannot be less than zero.")
                        return None
                    if not int_start < int_end:
                        print(f"{Fore.RED}[Mir-Us]   Wrong coordinates; start cannot be lower than end")
                        return None
                    for res in result:
                        for key in res.genome_coordinates:
                            for coord in res.genome_coordinates[key]:
                                try:
                                    f_start = int(coord[0])
                                    f_stop = int(coord[1])
                                except:
                                    continue
                                if (int_start <= f_start < int_end) and (int_start < f_stop <= int_end) \
                                        and not utils._exists(temp_result, res):
                                    temp_result.append(res)
                    result = temp_result
            if not end and start:
                temp_result = []
                if first:
                    int_start = int(start)
                    if int_start < 0:
                        print(f"{Fore.RED}[Mir-Us]   Incorrect 'start' or 'end' value; 'start' or 'end' cannot be less than zero.")
                        return None
                    for mi in self._miRNAs_ID:
                        for key in self._miRNAs_ID[mi].genome_coordinates:
                            for coord in self._miRNAs_ID[mi].genome_coordinates[key]:
                                try:
                                    f_start = int(coord[0])
                                except:
                                    continue
                                if (int_start <= f_start) and not utils._exists(
                                        result, self._miRNAs_ID[mi]):
                                    result.append(self._miRNAs_ID[mi])
                    first = False
                elif not first:
                    int_start = int(start)
                    if int_start < 0:
                        print(f"{Fore.RED}[Mir-Us]   Incorrect 'start' or 'end' value; 'start' or 'end' cannot be less than zero.")
                        return None
                    # count = 0
                    for res in result:
                        # count += 1
                        for key in res.genome_coordinates:
                            for coord in res.genome_coordinates[key]:
                                try:
                                    f_start = int(coord[0])
                                    # print(f"Current {(count)}: {f_start}")
                                except:
                                    continue
                                if int_start <= f_start and not utils._exists(temp_result, res):
                                    temp_result.append(res)
                                    # print(f"Kept:{f_start}")
                    result = temp_result
                # print("Only start")
            if not start and end:
                temp_result = []
                if first:
                    int_end = int(end)
                    if int_end < 0:
                        print(f"{Fore.RED}[Mir-Us]   Incorrect 'start' or 'end' value; 'start' or 'end' cannot be less than zero.")
                        return None
                    for mi in self._miRNAs_ID:
                        for key in self._miRNAs_ID[mi].genome_coordinates:
                            for coord in self._miRNAs_ID[mi].genome_coordinates[key]:
                                try:
                                    f_stop = int(coord[1])
                                except:
                                    continue
                                if (int_end >= f_stop) and not utils._exists(
                                        result, self._miRNAs_ID[mi]):
                                    result.append(self._miRNAs_ID[mi])
                    first = False
                elif not first:
                    int_end = int(end)
                    if int_end < 0:
                        print(f"{Fore.RED}[Mir-Us]   Incorrect 'start' or 'end' value; 'start' or 'end' cannot be less than zero.")
                        return None
                    # count = 0
                    for res in result:
                        # count += 1
                        for key in res.genome_coordinates:
                            for coord in res.genome_coordinates[key]:
                                try:
                                    f_stop = int(coord[1])
                                    # print(f"Current {(count)}: {f_start}")
                                except:
                                    continue
                                if int_end >= f_stop and not utils._exists(temp_result, res):
                                    temp_result.append(res)
                                    # print(f"Kept:{f_start}")
                    result = temp_result
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
                        if not utils._exists(result, mi):
                            result.append(self._miRNAs_ID[mi])
                except:
                    continue
            dict_result["prec_id-search"] = result
        if name:
            result = []
            try:
                if not utils._exists(result, self._miRNAs_ID[self._matures_name[name]]):
                    result.append(self._miRNAs_ID[self._matures_name[name]])
            except:
                pass
            dict_result["name-search"] = result
        if not start and not end and not chr and not strand and organism_name:
            result = []
            for mi in self._miRNAs_ID:
                if (self._miRNAs_ID[mi].organism == organism_name) and not utils._exists(result, self._miRNAs_ID[mi]):
                    result.append(self._miRNAs_ID[mi])
            dict_result["organism-search"] = result
        if tax_level:
            result = []
            for prec in self._precursors_ID:
                if tax_level in self._precursors_ID[prec].taxonomy:
                    mi_id = self._precursors_ID[prec].miRNAs
                    for elem in mi_id:
                        if not utils._exists(result, self._miRNAs_ID[elem]):
                            result.append(self._miRNAs_ID[elem])
            dict_result["tax-search"] = result
        if not start and not end and not organism_name and not strand and chr:
            result = []
            for mi in self._miRNAs_ID:
                if (chr in self._miRNAs_ID[mi].chromosome) and not utils._exists(result, self._miRNAs_ID[mi]):
                    result.append(self._miRNAs_ID[mi])
            dict_result["chr-search"] = result
        if not start and not end and not organism_name and not chr and strand:
            if not utils._check_strand(strand):
                print(f"{Fore.RED}[Mir-Us]   Incorrect 'strand' value; 'strand' can only be '+' or '-'")
                return None
            result = []
            for mi in self._miRNAs_ID:
                if (strand in self._miRNAs_ID[mi].strand) and not utils._exists(result, self._miRNAs_ID[mi]):
                    result.append(self._miRNAs_ID[mi])
            dict_result["strand-search"] = result
        if len(dict_result) >= 1 and not bool([res for res in dict_result.values() if res != []]):
            return None
        if not len(dict_result) > 1:
            # print(f"print: {dict_result}")
            # print("less than 1 key")
            # if list(dict_result.values())[0]:
            #if bool(dict_result):
            if dict_result.values():
                # print(dict_result.values())
                # print("if values")
                output = list(dict_result.values())[0]
                return output, len(output)
            else:
                return None
        # if len(dict_result) >= 1 and not bool([res for res in dict_result.values() if res != []]):
        #     return None
        print(f"{Fore.YELLOW}[Mir-Us]   Some of the given criteria were contradicting for the search system. Because of"
              f" that, the results are returned in a dictionary, where contradicting results are separated into"
              f" different search types. This search consist of (keys of generated dictionary): ")
        for key in dict_result.keys():
            print(f"{Fore.BLUE}'{key}': {len(dict_result[key])} results", sep="\n")
        # print([key for key in dict_result.keys()], sep="\n")
        return dict_result, sum([len(res) for res in list(dict_result.values())])

    @utils.time_this
    def find_cluster(self, mirna_id: str = None, prec_id: str = None, search_type="up-downstream", range=None, verbose=False):
        """Returns all miRNAs present within given range from given miRNA in affiliated organism genome.

        Args:
            mirna_id (str): miRNA ID
            prec_id (str): Precursor ID
            search_type (Union[str, int]): Search type which determines how to search genomic space ([**additional information**](#additional-information)).
            range (str): Length of genome to be searched
            verbose (bool): A flag, which allows or disallows showing search details (number of returned elements, time of execution, errors)

        Returns:
            Optional[list[MiRNA]]: miRNA objects matching given criteria or `None` if no results are found

        ## Additional information:
        ???+ info "Cluster search types"

            | Search type            | Explanation                                             |
            | ---------------------- | ------------------------------------------------------- |
            | `up-downstream` or `0` | The search is conducted below and above given miRNA position (towards both ends) within given range. The search is inclusive on both ends of the specified range.|
            | `upstream` or `1`      | The search is conducted only below given miRNA position (towards 5' end) within given range. The search is inclusive on both ends of the specified range.|
            | `downstream` or `2`    | The search is conducted only above given miRNA position (towards 3' end) within given range. The search is inclusive on both ends of the specified range.|
        """

        result = []

        def search_prec2(self, start, org, range):
            int_start = 0
            int_end = 0
            count = 0
            if search_type == "up-downstream" or search_type == 0:
                int_start = start - range
                # print(int_start)
                int_end = start + range
                # print(int_end)
            if search_type == "downstream" or search_type == 2:
                int_start = start
                int_end = start + range
            if search_type == "upstream" or search_type == 1:
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
                                int_start < f_stop <= int_end) and not utils._exists(result, self._precursors_ID[prec]):
                            count += 1
                            print(f"{count}, {self._precursors_ID[prec].ID}: {coord}")
                            result.append(self._precursors_ID[prec])

        def search_mirna2(self, start, org, range):
            int_start = 0
            int_end = 0
            count = 0
            if search_type == "up-downstream" or search_type == 0:
                int_start = start - range
                print(f"new_start: {int_start}")
                int_end = start + range
                print(f"new_end: {int_end}")
            elif search_type == "upstream" or search_type == 1:
                int_start = start
                int_end = start + range
            elif search_type == "downstream" or search_type == 2:
                int_start = start - range
                # print(int_start)
                int_end = start
            else:
                print(f"{Fore.RED}[Mir-Us]   Incorrect search type; possible search types are:\n "
                      f" - 'up-downstream' or '0'\n  - 'upstream' or '1'\n  - 'downstream' or '2'")
                return None
            for prec in self._precursors_ID:
                if self._precursors_ID[prec].organism == org:
                #if self._miRNAs_ID[mirna_id].organism == org and mirna_id in self._precursors_ID[prec].miRNAs:
                    for coord in self._precursors_ID[prec].genome_coordinates:
                        try:
                            # print(coord)
                            f_start = int(coord[0])
                            f_stop = int(coord[1])
                        except:
                            continue
                        if (int_start <= f_start < int_end) and (
                                int_start < f_stop <= int_end) and not utils._exists(result, self._precursors_ID[prec]):
                            count += 1
                            print(f"{count}, {self._precursors_ID[prec].ID}: {coord}")
                            result.append(self._precursors_ID[prec])
        # print(f"Mirna: {mirna_id} Prec: {prec_id}")
        if (mirna_id is not None or "") and (prec_id is not None or ""):
            print(f"{Fore.RED}[Mir-Us]   Contradicting actions; clusters cannot be searched between different types of "
                  f"objects.\nPlease, search MiRNA and Precursor objects separately.")
            return None
        elif mirna_id:
            if int(range) < 0:
                print(f"{Fore.RED}[Mir-Us]   Incorrect 'range' value; range cannot be less than zero.")
                return None
            org = self._miRNAs_ID[mirna_id].organism
            #print(org)
            # for key in self._miRNAs_ID[mirna_id].genome_coordinates:
            #     for coord in self._miRNAs_ID[mirna_id].genome_coordinates[key]:
            for coord in self._precursors_ID[self._miRNAs_ID[mirna_id].precursors[0]].genome_coordinates:
                start_position = int(coord[0])
                # print(f"org_start: {start_position}")
                search_mirna2(self, start_position, org, int(range))
        elif prec_id:
            if int(range) < 0:
                print(f"{Fore.RED}[Mir-Us]   Incorrect 'range' value; range cannot be less than zero.")
                return None
            org = self._precursors_ID[prec_id].organism
            for coord in self._precursors_ID[prec_id].genome_coordinates:
                start_position = int(coord[0])
                search_prec2(self, start_position, org, int(range))
        if not result:
            return None

        return result, len(result)

    @utils.time_this
    def get_tree(self, tax_path=None, verbose=False):
        """Returns taxonomy tree, where each taxonomy level is a dictionary (nested dictionaries as access to another
        taxonomy level and list of organisms at particular taxonomy level if has any)

        Args:
            tax_path (list[str]): Taxonomy path which indicates tree accession (tree is sliced from positon this path points to)
            verbose (bool): A flag, which allows or disallows showing search details (number of returned elements, time of execution, errors)

        Returns:
            Optional[dict]: Dictionary structured as taxonomy tree of organisms present in the base
            (equivalent to [Browse miRBase by species](https://www.mirbase.org/cgi-bin/browse.pl){: target="_blank"}) or its fragment.
            Returns `None` if wrong taxonomy path is given.
        """

        tax_dct = self._make_tax_dict()

        # structure def
        def tree():
            return dd(tree)

        # function to create tree (nested default dicts)
        def add_tree(t, keys):
            for key in keys:
                t = t[key]

        # function to access keys
        def access_tree(t, keys):
            for key in keys:
                t = t[key]
            return t

        # create tree with taxonomy
        org_tree = tree()
        for key in tax_dct:
            add_tree(org_tree, tax_dct[key])

        # add organisms
        for key in tax_dct:
            access_tree(org_tree, tax_dct[key])["!organism"] = [key]
        for key in tax_dct:
            if key not in access_tree(org_tree, tax_dct[key])["!organism"]:
                access_tree(org_tree, tax_dct[key])["!organism"].append(key)

        # sort in organism lists
        for key in tax_dct:
            access_tree(org_tree, tax_dct[key])["!organism"].sort()

        if tax_path is not None:
            try:
                tree_slice = dict(access_tree(org_tree, tax_path))
                if not tree_slice:
                    return None
                #print(json.dumps(tree_slice, indent=4, sort_keys=True))
                return tree_slice, 1
            except:
                return None
        else:
            #print(json.dumps(dict(org_tree), indent=4, sort_keys=True))
            return dict(org_tree), 1

    @utils.time_this
    def high_conf(self, mirna_obj=None, prec_obj=None, verbose=False):
        """From the given MiRNA or Precursor objects returns only those that are of a high confidence.

        Args:
            mirna_obj (list[MiRNA]): List of MiRNA objects, that can be retrieved using `get_mirna()` function
            prec_obj (list[Precursor]): List of Precursor objects, that can be retrieved using `get_precursor()` function
            verbose (bool): A flag, which allows or disallows showing search details (number of returned elements, time of execution, errors)

        Returns:
            Optional[list]: Filtered MiRNA or Precursor objects that are of a high confidence. Returns `None` if wrong
            objects are given or if none of the given objects are of a high confidence.
        """
        # if isinstance(mirna_obj, miObject.MiRNA):
        #     mirna_obj = [mirna_obj]
        # if isinstance(prec_obj, miObject.Precursor):
        #     prec_obj = [prec_obj]
        result = []
        if prec_obj:
            try:
                result = [prec for prec in prec_obj if prec.high_confidence is True]
            except:
                pass
        if mirna_obj:
            try:
                for mirna in mirna_obj:
                    for prec in mirna.precursors:
                        if self._precursors_ID[prec].high_confidence is True:
                            result.append(mirna)
            except:
                pass
        if not result:
            return None
        return list(set(result)), len(result)

    def dump_sequences(self, mirna_obj=None, prec_obj=None, filepath="", verbose=False):
        """Function writes sequences from given MiRNA or Precursor objects to a file in a FASTA format.

        Args:
            mirna_obj (list[MiRNA]): List of MiRNA objects, that can be retrieved using `get_mirna()` function
            prec_obj (list[Precursor]): List of Precursor objects, that can be retrieved using `get_precursor()` function
            filepath (str): String with the name of the file to which sequences will be saved.
            verbose (bool): A flag, which allows or disallows showing search details (number of returned elements, time of execution, errors)
        """
        start = timer()
        records = []
        if prec_obj:
            try:
                records = [SeqRecord(Seq(prec.precursor_sequence), id=prec.ID, name=prec.name, description=f"{prec.name}") for prec in prec_obj]
                # with open(f"{filepath}.fasta", "w") as handle:
                #     SeqIO.write(records, handle, "fasta")
                # for prec in prec_obj:
                #     record = SeqRecord(
                #         Seq(prec.precursor_sequence),
                #         id=prec.ID,
                #         name=prec.name
                #     )
            except:
                pass
        if mirna_obj:
            pass
        try:
            with open(f"{filepath}.fasta", "w") as handle:
                SeqIO.write(records, handle, "fasta")
        except IOError:
            print(f"{Fore.RED}[Mir-Us]   'dump_sequences' was unsuccessful.")
            return
        end = timer()
        runtime = end - start
        if verbose:
            print(f"{Fore.GREEN}[Mir-Us]   'dump_sequences' wrote {len(records)} sequences in {runtime:.6f} seconds")


class MiRLoad(MiRBase):
    """
    Parsing/loading class, which contains:
    - all functions to parse information directly from MiRBase using ftp
    - only newest version of MiRBase is supported
    """

    def __init__(self, version):
        super().__init__(version)

    def load_organisms(self, file_path):
        """
        Parses file with organisms (for 'current' database version - organisms.txt.gz)

        :param file_path: ftp path to files
        :type file_path: str
        :return: Creates namedtuples of organisms
        """

        def parse_organism_file(organisms_file):
            Entrez.email = "Your.Name.Here@example.org"
            for line in organisms_file:
                if self._miRBase_version is not "CURRENT" and int(self._miRBase_version) < 19:
                    line = line.decode('utf-8')
                if not line.startswith("#"):
                    tmp = line.split('\t')
                    self._org_sh[tmp[0]] = tmp[2]
                    if self._miRBase_version is not "CURRENT" and int(self._miRBase_version) < 20:
                        ent_handle = Entrez.esearch(db="taxonomy", retmax=10, term=tmp[2])
                        ent_record = Entrez.read(ent_handle)
                        try:
                            org = self._Organism(tmp[0], tmp[1], tmp[2], tmp[3].rstrip(), str(ent_record["IdList"][0]))
                        except:
                            org = self._Organism(tmp[0], tmp[1], tmp[2], tmp[3].rstrip(), None)
                        # print(f"{getattr(org, 'name')}: {getattr(org, 'taxid')}")
                        self._organisms.append(org)
                    else:
                        org = self._Organism(tmp[0], tmp[1], tmp[2], tmp[3], tmp[4].rstrip())
                        self._organisms.append(org)

        with urllib.request.urlopen(file_path) as organisms_file_url:
            if self._miRBase_version is not "CURRENT" and int(self._miRBase_version) < 19:
                parse_organism_file(organisms_file_url)
            else:
                with gzip.open(organisms_file_url, mode='rt') as organisms_file:
                    parse_organism_file(organisms_file)

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
                                self._matures_name[miRNA_vals[3]] = miRNA_vals[2]
                                # self._miRNAs_ID[miRNA_vals[2]].mature_sequence.append(self._miRNAs_ID[miRNA_vals[2]].get_mature_seq(seq_full, (miRNA_vals[0], miRNA_vals[1])))
                                # self._miRNAs_ID[miRNA_vals[2]] = mature.to_json()
                                # self._miRNAs_ID[miRNA_vals[2]] = mature.dumper()
                            else:
                                self._miRNAs_ID[miRNA_vals[2]].precursors.append(id_p)
                                if miRNA_vals[3] not in self._miRNAs_ID[miRNA_vals[2]].name:
                                    self._miRNAs_ID[miRNA_vals[2]].name.append(miRNA_vals[3])
                                self._miRNAs_ID[miRNA_vals[2]].mature_positions.append((miRNA_vals[0], miRNA_vals[1]))
                                self._miRNAs_ID[miRNA_vals[2]].evidence.append(miRNA_vals[4])
                                self._miRNAs_ID[miRNA_vals[2]].experiment.append(miRNA_vals[5])
                                self._miRNAs_ID[miRNA_vals[2]].end.append(e)
                                # self._miRNAs_ID[miRNA_vals[2]].mature_sequence.append(
                                #     self._miRNAs_ID[miRNA_vals[2]].get_mature_seq(seq_full,
                                #                                                   (miRNA_vals[0], miRNA_vals[1])))
                                if len(ref) > 1:
                                    pass
                                    # self._miRNAs_ID[miRNA_vals[2]].references.append(ref)
                                    # self._miRNAs_ID[miRNA_vals[2]].references.append(''.join(map(str, ref)))
                                else:
                                    if ''.join(map(str, ref)) not in self._miRNAs_ID[miRNA_vals[2]].references:
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

        iupac = ['a', 'c', 'g', 't', 'u', 'r', 'y', 's', 'w', 'k', 'm', 'b', 'd', 'h', 'v', 'n']
        start = 0
        end = 7
        x = 0
        with urllib.request.urlopen(file_path) as struct_file_gz:
            struct_file = gzip.open(struct_file_gz, mode='rt')
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
                    acgu_dt = ''.join(acgu_dt).lower().replace("-", "")
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
        # print(tax_dct)

        for prec in self._precursors_ID:
            p_id = self._precursors_ID[prec].ID
            p_org = self._precursors_ID[prec].organism
            self._organisms_of_prec[p_org].append(p_id)
            tax_name = tax_dct[p_org].split(';')[:-1]
            self._precursors_ID[prec].taxonomy = tax_name
            for tax_level in tax_name:
                if not (p_id in self._taxonomy_of_prec[tax_level]):
                    self._taxonomy_of_prec[tax_level].append(p_id)

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
                    sim_Alias = split_info_miRNA[1].split('=')[1]

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
                            self._miRNAs_ID[sim_Alias].chromosome.append(chr_miRNA)
                            self._miRNAs_ID[sim_Alias].strand.append(strand_seq)
                            self._miRNAs_ID[sim_Alias].genome_coordinates[derivative].append((start_seq, end_seq))
                        except:
                            continue
