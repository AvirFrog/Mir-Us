import miBase

#
# example usage of Mir-Us
# test examples of usage below
# some interesting outputs below test examples
# whole output is dumped into 'mirek_out.txt'

# initialising
m = miBase.MiRBase(version="CURRENT")

# actual test examples
print("----- get_organisms_list -----")
m.get_organisms_list()  # nothing will appear in output - just returning list of all organisms
print(m.get_organisms_list())  # prints list of all organisms

print("----- get_organisms_short -----")
m.get_organisms_short()  # nothing will appear in output - just returning list of all organism abbreviations
print(m.get_organisms_short())  # prints list of all organism abbreviations

print("----- get_organism-----")
m.get_organism()  # no results message
m.get_organism("Viruses")  # successful search message with elapsed time
m.get_organism("Dogs")  # no results message
print(m.get_organism())  # no results message and then printed 'None'
print(m.get_organism("Viruses"))  # successful search message with elapsed time and then printed list with results
print(m.get_organism("Dogs"))  # no results message and then printed 'None'

print("----- get_taxid-----")
m.get_taxid()  # no results message
m.get_taxid(["Amphimedon queenslandica", "Homo sapiens", "Chrysemys picta"])  # successful search message with elapsed time
m.get_taxid(["Homo", "Chrysemys picta"])  # successful search message with elapsed time
m.get_taxid(["Sapiens"])  # no results message
m.get_taxid("Chrysemys picta")  # successful search message with elapsed time
print(m.get_taxid())  # no results message and then printed 'None'
print(m.get_taxid(["Amphimedon queenslandica", "Homo sapiens", "Chrysemys picta"]))  # successful search message with elapsed time and printed dict with results
print(m.get_taxid(["Homo", "Chrysemys picta"]))  # successful search message with elapsed time and printed dict with results
print(m.get_taxid(["Sapiens"]))  # no results message and then printed 'None'
print(m.get_taxid("Chrysemys picta"))  # successful search message with elapsed time and printed dict with results

print("----- get_tax_level-----")
m.get_tax_level()  # no results message
m.get_tax_level(['Homo', "Amphimedon queenslandica"])  # successful search message with elapsed time
m.get_tax_level(['Homo'])  # no results message
m.get_tax_level("Homo sapiens")  # successful search message with elapsed time
print(m.get_tax_level())  # no results message and then printed 'None'
print(m.get_tax_level(['Homo', "Amphimedon queenslandica"]))  # successful search message with elapsed time and printed dict with results
print(m.get_tax_level(['Homo']))  # no results message and then printed 'None'
print(m.get_tax_level("Homo sapiens"))  # successful search message with elapsed time and printed dict with results

print("----- get_structure-----")
m.get_structure()  # no results message
m.get_structure(["MI0000001", "MI0016085"])  # successful search message with elapsed time
m.get_structure(["MI123456789"])  # no results message
m.get_structure("MI0000001")  # successful search message with elapsed time
print(m.get_structure())  # no results message and then printed 'None'
print(m.get_structure(["MI0000001", "MI0016085"]))  # successful search message with elapsed time and printed dict with results
print(m.get_structure(["MI123456789"]))  # no results message and then printed 'None'
print(m.get_structure("MI0000001"))  # successful search message with elapsed time and printed dict with results

print("----- get_references-----")
m.get_references()  # no results message
m.get_references(mirna_id=['MIMAT0000001', "MIMAT123"])  # successful search message with elapsed time
m.get_references(['MIMAT0000001', "MIMAT123"])  # successful search message with elapsed time
m.get_references(['MIMAT123'])  # no results message
m.get_references('MIMAT0000001')  # successful search message with elapsed time
m.get_references("gga-miR-7478-3p")  # no results message - mirna_id is default
m.get_references(mirna_name="gga-miR-7478-3p")  # successful search message with elapsed time
m.get_references(mirna_name=["gga-miR-7478-3p"])  # successful search message with elapsed time
m.get_references(["gga-miR-7478-3p"])  # no results message - mirna_id is default
m.get_references(['MIMAT0000001', "gga-miR-7478-3p"])  # successful search message with elapsed time; will only search through mirna_id
print(m.get_references())  # no results message and then printed 'None'
print(m.get_references(mirna_id=['MIMAT0000001', "MIMAT123"]))  # successful search message with elapsed time and printed list with results
print(m.get_references(['MIMAT0000001', "MIMAT123"]))  # successful search message with elapsed time and printed list with results
print(m.get_references(['MIMAT123']))  # no results message and then printed 'None'
print(m.get_references('MIMAT0000001'))  # successful search message with elapsed time and printed list with results
print(m.get_references("gga-miR-7478-3p"))  # no results message and then printed 'None'
print(m.get_references(mirna_name="gga-miR-7478-3p"))  # successful search message with elapsed time and printed list with results
print(m.get_references(mirna_name=["gga-miR-7478-3p"]))  # successful search message with elapsed time and printed list with results
print(m.get_references(["gga-miR-7478-3p"]))  # no results message and then printed 'None'
print(m.get_references(['MIMAT0000001', "gga-miR-7478-3p"]))  # successful search message with elapsed time and printed list with results

print("----- get_precursor-----")
m.get_precursor()  # no results message
m.get_precursor("MI0000001")  # successful search message with elapsed time
m.get_precursor(name="mmu-mir-21a")  # successful search message with elapsed time
m.get_precursor(id=["MI0000001", "MI0017717", "MI0000021"])  # successful search message with elapsed time
m.get_precursor(id=["MI9999000"])  # no results message
m.get_precursor("MI9999000")  # no results message
m.get_precursor("Cnidaria")  # no results message
m.get_precursor(start=186, end=11046)  # successful search message with elapsed time
m.get_precursor(start='186', end='11046')  # successful search message with elapsed time
m.get_precursor(strand='+')  # successful search message with elapsed time
m.get_precursor(organism_name="Homo sapiens")  # successful search message with elapsed time
m.get_precursor(mirna_id="MIMAT0000001")  # successful search message with elapsed time
m.get_precursor(tax_level="Cnidaria")  # successful search message with elapsed time
m.get_precursor(id=["MI9999000"], organism_name="Homo sapiens", chr='chrX')  # successful search message with elapsed time
print(m.get_precursor())  # no results message and then printed 'None
print(m.get_precursor("MI0000001"))  # successful search message with elapsed time and then pretty printed records
print(m.get_precursor(name="mmu-mir-21a"))  # successful search message with elapsed time and then pretty printed records
print(m.get_precursor(id=["MI0000001", "MI0017717", "MI0000021"]))  # successful search message with elapsed time and then pretty printed records
print(m.get_precursor(id=["MI9999000"]))  # no results message and then printed 'None
print(m.get_precursor("MI9999000"))  # no results message and then printed 'None
print(m.get_precursor("Cnidaria"))  # no results message and then printed 'None
print(m.get_precursor(start=186, end=11046))  # successful search message with elapsed time and then pretty printed records
print(m.get_precursor(start='186', end='11046'))  # successful search message with elapsed time and then pretty printed records
print(m.get_precursor(strand='+'))  # successful search message with elapsed time and then pretty printed records
print(m.get_precursor(organism_name="Homo sapiens"))  # successful search message with elapsed time and then pretty printed records
print(m.get_precursor(mirna_id="MIMAT0000001"))  # successful search message with elapsed time and then pretty printed records
print(m.get_precursor(tax_level="Cnidaria"))  # successful search message with elapsed time and then pretty printed records
print(m.get_precursor(id=["MI9999000"], organism_name="Homo sapiens", chr='chrX'))  # successful search message with elapsed time and then pretty printed records

print("----- get_mirna-----")
m.get_mirna()  # no results message
m.get_mirna(id=['MIMAT0000001', 'MIMAT123'])  # successful search message with elapsed time
m.get_mirna("MIMAT0000001")  # successful search message with elapsed time
m.get_mirna('MIMAT123')  # no results message
m.get_mirna("Hominidae")  # no results message
m.get_mirna(tax_level="Homo sapiens", start=11046, end=186)  # wrong coordinates message and then no results message
m.get_mirna(tax_level="Hominidae")  # successful search message with elapsed time
m.get_mirna(id=["MIMAT0000001"], organism_name="Brassica napus", chr='chr5')  # successful search message with elapsed time
print(m.get_mirna())  # no results message
print(m.get_mirna(id=['MIMAT0000001', 'MIMAT123']))  # successful search message with elapsed time
print(m.get_mirna("MIMAT0000001"))  # successful search message with elapsed time
print(m.get_mirna('MIMAT123'))  # no results message
print(m.get_mirna("Hominidae"))  # no results message
print(m.get_mirna(tax_level="Homo sapiens", start=11046, end=186))  # wrong coordinates message and then no results message
print(m.get_mirna(tax_level="Hominidae"))  # successful search message with elapsed time
print(m.get_mirna(id=["MIMAT0000001"], organism_name="Brassica napus", chr='chr5'))  # successful search message with elapsed time


# example results:
# instresting one ;)
#         Mature ID: MIMAT0029111
#         Mature name: gga-miR-7478-3p
#         Derivative: MI0024151
#         Mature sequence: cuacaugcacgggcaggugag
#         Mature positions: ('58', '78')
#         Evidence: experimental
#         Experiment: Illumina [1]
#         End: 3p
#         Chromosome: AADN04000787.1, AADN04001013.1, AADN04001058.1, AADN04001184.1, AADN04001195.1, AADN04001464.1, AADN04001515.1, AADN04001825.1, AADN04001939.1, AADN04002188.1, AADN04002235.1, AADN04002674.1, AADN04002678.1, AADN04003570.1, AADN04007859.1, AADN04009770.1, AADN04012884.1, AADN04016299.1, KQ759544.1
#         Genome coordinates: {'MI0024151': [('986', '1006'),
#                ('28464', '28484'),
#                ('20409', '20429'),
#                ('24105', '24125'),
#                ('17106', '17126'),
#                ('9726', '9746'),
#                ('12927', '12947'),
#                ('10860', '10880'),
#                ('14590', '14610'),
#                ('4521', '4541'),
#                ('14892', '14912'),
#                ('12904', '12924'),
#                ('5898', '5918'),
#                ('12564', '12584'),
#                ('3531', '3551'),
#                ('3948', '3968'),
#                ('6117', '6137'),
#                ('4875', '4895'),
#                ('1097', '1117')]}
#         Strand: -, -, -, -, +, +, +, -, -, -, -, +, -, +, +, +, +, -, -
#         References: 23034410

# BIG one
#         Mature ID: MIMAT0029118
#         Mature name: gga-miR-7482-5p, gga-miR-7482-5p, gga-miR-7482-5p, gga-miR-7482-5p, gga-miR-7482-5p, gga-miR-7482-5p
#         Derivative: MI0024162, MI0031115, MI0031116, MI0031117, MI0031118, MI0031119
#         Mature sequence: ccugggcuuguucacucaccagaga, ccugggcuuguucacucaccagaga, ccugggcuuguucacucaccagaga, ccugggcuuguucacucaccagaga, ccugggcuuguucacucaccagaga, ccugggcuuguucacucaccagaga
#         Mature positions: ('1', '25'), ('1', '25'), ('1', '25'), ('1', '25'), ('1', '25'), ('1', '25')
#         Evidence: experimental, experimental, experimental, experimental, experimental, experimental
#         Experiment: Illumina [1], Illumina [1], Illumina [1], Illumina [1], Illumina [1], Illumina [1]
#         End: 5p, 5p, 5p, 5p, 5p, 5p
#         Chromosome: chrZ, chrZ, chrZ, chrZ, chrZ, KQ759519.1
#         Genome coordinates: {'MI0024162': [('34300480', '34300504')],
#  'MI0031115': [('34315561', '34315585')],
#  'MI0031116': [('34330947', '34330971')],
#  'MI0031117': [('34395460', '34395484')],
#  'MI0031118': [('34410494', '34410518')],
#  'MI0031119': [('10680', '10704')]}
#         Strand: +, +, +, +, +, -
#         References: 23034410, 23034410, 23034410, 23034410, 23034410, 23034410, 23034410, 23034410, 23034410, 23034410, 23034410