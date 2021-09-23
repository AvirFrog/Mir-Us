import miBase
#from miBase import MiRBase._compile_indexes()
# import numpy as np
# from timeit import default_timer as timer

# example usage of Mir-Us
# test examples of usage below
# some interesting outputs below test examples
# whole output is dumped into 'mirek_out.txt'

# initialising
m = miBase.MiRBase(version="CURRENT")
#m._compile_indexes()

# l = list()
# l2 = list()
# lnp = np.asarray(l2)
#
# start = timer()
# for i in range(100000):
#     l.append(i)
# end = timer()
# print(f"Normal list: {end - start}")
#
# start = timer()
# for i in range(100000):
#     np.append(lnp, i)
# end = timer()
# print(f"Numpy array: {end - start}")
#
# start = timer()
# for i in l:
#     j = i + i
# end = timer()
# print(f"Normal list add: {end - start}")
#
# start = timer()
# for i in range(lnp.shape[0]):
#     j = lnp[i] + lnp[i]
# end = timer()
# print(f"Numpy array add: {end - start}")

#print(m.get_precursor(name="cel-let-7"))
#m.get_precursor(name=["cin-mir-4058", "cel-mir-4931", "aly-MIR395h", "dma-mir-423", "hsa-mir-3613"])

# # actual test examples
# print("----- get_organisms_list -----")
# organisms = m.get_organisms_list()  # nothing will appear in output - just returning list of all organisms
# print(m.get_organisms_list())  # prints list of all organisms
# gallus = [nt for nt in organisms if getattr(nt, "name") == "Gallus gallus"][0]
# print(gallus)
# print(getattr(gallus, "organism"))
# print(getattr(gallus, "name"))
# print(getattr(gallus, "tree"))
# print(getattr(gallus, "taxid"))
# shorts = m.get_organisms_short()
# gallus = shorts['gga']
# print(gallus)
#
# short_gal = list(shorts.keys())[list(shorts.values()).index("Gallus gallus")]
# print(short_gal)
#
# print("----- get_organisms_short -----")
# m.get_organisms_short()  # nothing will appear in output - just returning list of all organism abbreviations
# #print(m.get_organisms_short())  # prints list of all organism abbreviations
#
# print("----- get_organism-----")
# m.get_organism()  # no results message
# m.get_organism("Viruses")  # successful search message with elapsed time
# m.get_organism("Dogs")  # no results message
# print(m.get_organism())  # no results message and then printed 'None'
# print(m.get_organism("Viruses"))  # successful search message with elapsed time and then printed list with results
# print(m.get_organism("Dogs"))  # no results message and then printed 'None'
print(m.get_organism("Aves"))

# print("----- get_taxid-----")
# m.get_taxid()  # no results message
# m.get_taxid(["Amphimedon queenslandica", "Homo sapiens", "Chrysemys picta"])  # successful search message with elapsed time
# m.get_taxid(["Homo", "Chrysemys picta"])  # successful search message with elapsed time
# m.get_taxid(["Sapiens"])  # no results message
#m.get_taxid("Chrysemys picta")  # successful search message with elapsed time
# #print(m.get_taxid())  # no results message and then printed 'None'
# print(m.get_taxid(["Amphimedon queenslandica", "Homo sapiens", "Chrysemys picta"]))  # successful search message with elapsed time and printed dict with results
# print(m.get_taxid(["Homo", "Chrysemys picta"]))  # successful search message with elapsed time and printed dict with results
# #print(m.get_taxid(["Sapiens"]))  # no results message and then printed 'None'
#print(m.get_taxid("Chrysemys picta"))  # successful search message with elapsed time and printed dict with results
print(m.get_taxid(["Homo", "Chrysemys picta", "Amphimedon queenslandica"]))
single = m.get_taxid("Gallus gallus")
print(single)
print(single["Gallus gallus"])

#
# print("----- get_tax_level-----")
# m.get_tax_level()  # no results message
# m.get_tax_level(['Homo', "Amphimedon queenslandica"])  # successful search message with elapsed time
# m.get_tax_level(['Homo'])  # no results message
# m.get_tax_level("Homo sapiens")  # successful search message with elapsed time
# print(m.get_tax_level())  # no results message and then printed 'None'
print(m.get_tax_level(['Homo', "Amphimedon queenslandica"]))  # successful search message with elapsed time and printed dict with results
# print(m.get_tax_level(['Homo']))  # no results message and then printed 'None'
print(m.get_tax_level("Gallus gallus"))  # successful search message with elapsed time and printed dict with results

# print("----- get_structure-----")
# m.get_structure()  # no results message
# m.get_structure(["MI0000001", "MI0016085"])  # successful search message with elapsed time
# m.get_structure(["MI123456789"])  # no results message
# m.get_structure("MI0000001")  # successful search message with elapsed time
# m.get_structure(name=["mmu-mir-21a", "hsa-mir-3612"])  # successful search message with elapsed time
# print(m.get_structure())  # no results message and then printed 'None'
# print(m.get_structure(["MI0000001", "MI0016085"]))  # successful search message with elapsed time and printed dict with results
# print(m.get_structure(["MI123456789"]))  # no results message and then printed 'None'
#print(m.get_structure("MI0000001"))  # successful search message with elapsed time and printed dict with results
# print(m.get_structure(name=["mmu-mir-21a", "hsa-mir-3612"]))  # successful search message with elapsed time and printed dict with results

# print("----- get_references-----")
# m.get_references()  # no results message
# m.get_references(mirna_id=['MIMAT0000001', "MIMAT123"])  # successful search message with elapsed time
# m.get_references(['MIMAT0000001', "MIMAT123"])  # successful search message with elapsed time
# m.get_references(['MIMAT123'])  # no results message
# m.get_references('MIMAT0000001')  # successful search message with elapsed time
# m.get_references("gga-miR-7478-3p")  # no results message - mirna_id is default
# m.get_references(mirna_name="gga-miR-7478-3p")  # successful search message with elapsed time
# m.get_references(mirna_name=["gga-miR-7478-3p"])  # successful search message with elapsed time
# m.get_references(["gga-miR-7478-3p"])  # no results message - mirna_id is default
# m.get_references(['MIMAT0000001', "gga-miR-7478-3p"], link=True)  # successful search message with elapsed time; will only search through mirna_id
# m.get_references(prec_id=["MI0000001", "MI0017717", "MI0000021"])  # successful search message with elapsed time
# print(m.get_references())  # no results message and then printed 'None'
# print(m.get_references(mirna_id=['MIMAT0000001', "MIMAT123"]))  # successful search message with elapsed time and printed list with results
# print(m.get_references(['MIMAT0000001', "MIMAT123"]))  # successful search message with elapsed time and printed list with results
# print(m.get_references(['MIMAT123']))  # no results message and then printed 'None'
# print(m.get_references('MIMAT0000001'))  # successful search message with elapsed time and printed list with results
# print(m.get_references("gga-miR-7478-3p"))  # no results message and then printed 'None'
# print(m.get_references(mirna_name="gga-miR-7478-3p"))  # successful search message with elapsed time and printed list with results
# print(m.get_references(mirna_name=["gga-miR-7478-3p"]))  # successful search message with elapsed time and printed list with results
# print(m.get_references(["gga-miR-7478-3p"]))  # no results message and then printed 'None'
# print(m.get_references(['MIMAT0000001', "gga-miR-7478-3p"], link=True))  # successful search message with elapsed time and printed list with results
# print(m.get_references(prec_id=["MI0000001", "MI0017717", "MI0000021"]))  # successful search message with elapsed time and printed list with results

# print("----- get_precursor-----")
# m.get_precursor()  # no results message
# m.get_precursor("MI0000001")  # successful search message with elapsed time
# m.get_precursor(name="mmu-mir-21a")  # successful search message with elapsed time
# m.get_precursor(id=["MI0000001", "MI0017717", "MI0000021"])  # successful search message with elapsed time
# m.get_precursor(id=["MI9999000"])  # no results message
# m.get_precursor("MI9999000")  # no results message
# m.get_precursor("Cnidaria")  # no results message
# m.get_precursor(start=186, end=11046)  # successful search message with elapsed time
# m.get_precursor(start='186', end='11046')  # successful search message with elapsed time
# m.get_precursor(strand='+')  # successful search message with elapsed time
# m.get_precursor(organism_name="Homo sapiens")  # successful search message with elapsed time
# m.get_precursor(mirna_id="MIMAT0000001")  # successful search message with elapsed time
# m.get_precursor(tax_level="Cnidaria")  # successful search message with elapsed time
# m.get_precursor(id=["MI9999000"], organism_name="Homo sapiens", chr='chrX')  # successful search message with elapsed time
# #print(m.get_precursor())  # no results message and then printed 'None
# #print(m.get_precursor("MI0000001"))  # successful search message with elapsed time and then pretty printed records
# #print(m.get_precursor(name="mmu-mir-21a"))  # successful search message with elapsed time and then pretty printed records
# #print(m.get_precursor(id=["MI0000001", "MI0017717", "MI0000021"]))  # successful search message with elapsed time and then pretty printed records
# #print(m.get_precursor(id=["MI9999000"]))  # no results message and then printed 'None
# #print(m.get_precursor("MI9999000"))  # no results message and then printed 'None
# #print(m.get_precursor("Cnidaria"))  # no results message and then printed 'None
# #print(m.get_precursor(start=186, end=11046))  # successful search message with elapsed time and then pretty printed records
# #print(m.get_precursor(start='186', end='11046'))  # successful search message with elapsed time and then pretty printed records
# #print(m.get_precursor(strand='+'))  # successful search message with elapsed time and then pretty printed records
# #print(m.get_precursor(organism_name="Homo sapiens"))  # successful search message with elapsed time and then pretty printed records
# #print(m.get_precursor(mirna_id="MIMAT0000001"))  # successful search message with elapsed time and then pretty printed records
# #print(m.get_precursor(tax_level="Cnidaria"))  # successful search message with elapsed time and then pretty printed records
# #print(m.get_precursor(id=["MI9999000"], organism_name="Homo sapiens", chr='chrX'))  # successful search message with elapsed time and then pretty printed records
#
# print("----- get_mirna-----")
# m.get_mirna()  # no results message
# m.get_mirna(id=['MIMAT0000001', 'MIMAT123'])  # successful search message with elapsed time
# m.get_mirna("MIMAT0000001")  # successful search message with elapsed time
# m.get_mirna('MIMAT123')  # no results message
# m.get_mirna("Hominidae")  # no results message
# m.get_mirna(tax_level="Homo sapiens", start=11046, end=186)  # wrong coordinates message and then no results message
# m.get_mirna(tax_level="Hominidae")  # successful search message with elapsed time
# m.get_mirna(id=["MIMAT0000001"], organism_name="Brassica napus", chr='chr5')  # successful search message with elapsed time
# #print(m.get_mirna())  # no results message
# #print(m.get_mirna(id=['MIMAT0000001', 'MIMAT123']))  # successful search message with elapsed time
#print(m.get_mirna("MIMAT0000001"))  # successful search message with elapsed time
# #print(m.get_mirna('MIMAT123'))  # no results message
# #print(m.get_mirna("Hominidae"))  # no results message
# #print(m.get_mirna(tax_level="Homo sapiens", start=11046, end=186))  # wrong coordinates message and then no results message
# #print(m.get_mirna(tax_level="Hominidae"))  # successful search message with elapsed time
# #print(m.get_mirna(id=["MIMAT0000001"], organism_name="Brassica napus", chr='chr5'))  # successful search message with elapsed time

#print(m.get_precursor("MI0000001"))
#print(m.get_mirna(organism_name="Caenorhabditis elegans"))
#print(m.find_cluster(prec_id="MI0000001", range="1000000"))
#print(m.get_precursor(id=["MI0026845", "MI0026846", "MI0026849"]))
#print(m.find_cluster(mirna_id="MIMAT0032841", range="1000000"))

# to ten dziwny przypadek
# m.find_cluster(mirna_id="MIMAT0032841", range="1000000")
# m.find_cluster(mirna_id="MIMAT0032841", search_type="upstream", range="1000000")
# m.find_cluster(mirna_id="MIMAT0032841", search_type="downstream", range="1000000")
#
# m.find_cluster(mirna_id="MIMAT0032841", range="100000")
# m.find_cluster(mirna_id="MIMAT0032841", search_type="upstream", range="100000")
# m.find_cluster(mirna_id="MIMAT0032841", search_type="downstream", range="100000")
#
# m.find_cluster(prec_id="MI0000001", range="100000")
# m.find_cluster(prec_id="MI0000001", search_type="upstream", range="100000")
# m.find_cluster(prec_id="MI0000001", search_type="downstream", range="100000")
#
# m.find_cluster(mirna_id="MIMAT0001641", range="100000")
# m.find_cluster(mirna_id="MIMAT0001641", search_type="upstream", range="100000")
# m.find_cluster(mirna_id="MIMAT0001641", search_type="downstream", range="100000")

# m.find_cluster(mirna_id="MIMAT0026631", range="1000")
# m.find_cluster(mirna_id="MIMAT0026631", search_type="upstream", range="1000")
# m.find_cluster(mirna_id="MIMAT0026631", search_type="downstream", range="1000")

#m.find_cluster(prec_id="MI0000072", range="10000")

#print(m.get_mirna(id="MIMAT0032841"))
#print(m.get_mirna(id="MIMAT0029111"))
#m.get_mirna(organism_name="Homo sapiens", id="MIMAT0029111")
# m.get_mirna(prec_id=["MI0000223", "MI0000182", "MI0000060", "MI000000"])
#print(m.get_mirna(prec_id=["MI0000223", "MI0000182", "MI0000060", "MI000000"]))
# m.get_precursor(mirna_id=["MIMAT0050213", "MIMAT0050254", "MIMAT0050056", "MIMAT0049408"])
#print(m.get_precursor(mirna_id=["MIMAT0050213", "MIMAT0050254", "MIMAT0050056", "MIMAT0049408"]))

#print(m.get_mirna(prec_id=["MI0000071"]))
print(m.get_precursor(prec_id=["MI0000071"]))
#print(m.get_mirna(chr="LSRX01000794.1", organism_name="Symbiodinium microadriaticum", start="144000"))
# obj = m.get_mirna(chr="chrX", organism_name="Homo sapiens", start="153300000", mirna_id="MIMAT0050213", tax_level="Viruses")
#m.get_mirna(organism_name="Homo sapiens", mirna_id="MIMAT0050213")
#print(type(obj))
#print(obj)
# print(m._miRNAs_ID["MIMAT0009754"])
#print(m._precursors_ID["MI0017682"])
#m.get_tree()
#print(m.get_mirna("MIMAT0000002")[0].__dict__)
# print(m.get_mirna("MIMAT0000002"))
# print(m.get_mirna("MIMAT0050213"))
# print(m.get_mirna("MI9999999999999"))
# m.get_mirna(chr="chrX", organism_name="Homo sapiens", start="153300000", mirna_id="MI9999999999999", tax_level="Viruses")
# m.get_mirna(chr="gowno", organism_name="Homo sapiens", start="153300000", mirna_id="MI9999999999999", tax_level="Virusez")
# tree = m.get_tree()

# gallus_mirna = m.get_mirna(organism_name="Gallus gallus", chr="chr8", strand="-", start="6000000")
# print(gallus_mirna)
# sample_gallus = m.get_mirna(mirna_id=["MIMAT0001185", "MIMAT0025825", "MIMAT0007451"])
# print(sample_gallus)
# struct_gallus = m.get_structure(id=[precursor for mi_obj in sample_gallus for precursor in mi_obj.precursor])
# print(struct_gallus)
# ref_gallus = m.get_references(mirna_id=[id for mi_obj in sample_gallus for id in mi_obj.mature_ID], link=True)
# print(ref_gallus)
# print(type(ref_gallus))

cluster_gallus = m.find_cluster(prec_id="MI0007558", range="10000")
print(cluster_gallus)

print(m.find_cluster(mirna_id="MIMAT0031077", range="10000"))

# print(m._precursors_ID["MI0005093"].genome_coordinates)
# print(m.get_precursor(prec_id=["MI0000001", "MI0017717", "MI0000021"]))
# m.get_precursor(prec_id=["MI9999000"], organism_name="Homo sapiens", chr='chrX')
# print(m.get_precursor(name="mghv-mir-M1-2"))
# #tree = m.get_tree(["Deuterostoma"])
# #print(tree["!organism"])
# tree = m.get_tree(["Metazoa", "Bilateria", "Deuterostoma", "Chordata", "Vertebrata"])
# print(tree["Aves"]["!organism"])
# aves = m.get_tree(["Metazoa", "Bilateria", "Deuterostoma", "Chordata", "Vertebrata", "Aves"])["!organism"]
# print(aves)
# #gallus_prec = m.get_precursor(organism_name="Gallus gallus")
# gallus_prec = m.get_precursor(organism_name="Gallus gallus", chr="chr8", strand="-", start="6000000")
# print(gallus_prec)
# print(m.get_precursor('MI0005093'))
# sample_gallus = m.get_precursor(prec_id=["MI0007316", "MI0022537", "MI0007409"])
# print(sample_gallus)
# print(sample_gallus[0].chromosome)
# print(sample_gallus[0].genome_coordinates)

print(m.get_mirna(mirna_id="MIMAT0031077"))

#print(dict(tree["Metazoa"]))
#m.get_mirna(chr="II", organism_name="Caenorhabditis elegans")
#m.get_mirna(chr="X", organism_name="Homo sapiens")
# print(m.get_mirna(mirna_id="MIMAT0002614"))
# print(m.get_mirna(mirna_id="MIMAT0029118"))
# print(m.get_precursor(id="MI0005013"))
# print(m.get_precursor(id="MI0003024"))
# print(m.get_mirna(mirna_id=["MIMAT0015054", "MIMAT0015056"]))
# print(m.get_mirna(name="hsa-miR-3179"))
# print(m.get_precursor(id="MI0008605"))
# print(m.get_precursor(name="ptr-mir-30c-1"))
# print(len(m._precursors_name))
# print(len(m._matures_name))
# print(m.get_mirna(chr="chrX", organism_name="Homo sapiens", start="153300000"))
# print(m.get_precursor("MI0015972"))
# print(m.get_mirna("MIMAT0005829"))
# hconfs = 0
# for record in m._precursors_ID:
#     if m._precursors_ID[record].high_confidence is True:
#         hconfs += 1
# print(hconfs)  # 3320 for 22.1
#                 # 2162 for 22
#                 # 1996 for 21
#                 # 0 from 20 and lower
# print(m.get_precursor("MI0001370"))
# print(type(m._precursors_ID["MI0008605"].genome_coordinates[0]))
# print(m.get_organisms_list())
# print(m.get_tax_level("Homo sapiens"))
#m.get_tree()
#uguaaacauccuacacucucagc                            agauacuguaaacauccuacacucucagcuguggaaaguaagaaagcugggagaaggcuguuuacucuuucu
#uguaaacauccuacacucucagc
#uguaguguguguaaacauccuac


#print(m._precursors_name)
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
