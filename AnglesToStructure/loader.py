import os

from Bio import *

from Bio.PDB import MMCIFParser


"""
loads s CIF file

Boris Starkov(c)
"""

data_folder = '../Data'

file_path = os.path.join(data_folder,'Residue_samples_CIF', 'ala.cif')

parser = MMCIFParser()



structure = parser.get_structure('hey', file_path)

print(structure)
