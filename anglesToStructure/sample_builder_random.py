import os
import sys
import gzip

from Bio.PDB import *
from Bio import SeqIO
from Bio.SeqUtils import seq1

import warnings

"""

Makes a set of structure forms of aminoacids just by taking first

Boris Starkov(c)

"""
warnings.filterwarnings('ignore')

data_folder = '../Data'
PDB_files_directory = os.path.join(data_folder, 'PDB_files')

result_folder = os.path.join(data_folder, 'Residue_samples')
if not os.path.exists(result_folder):
	os.makedirs(result_folder)


parser = PDBParser()

writer = PDBIO()

used = set()

#go through all the files
for root, dirs, files in os.walk(PDB_files_directory):
	for file_name in files:

		file_path = os.path.join(root, file_name) # path to current file

		with gzip.open(file_path, 'rb') as file:
			structure = parser.get_structure('', file)

		#### checks if a structure fits


		## bad resolution
		if (structure.header['resolution'] == None or structure.header['resolution'] > 1): #bad resolution
			print('bad res')
			continue 

		## more than 1 model
		models_counter = 0
		for model in structure.get_models():
			models_counter += 1
			if (models_counter > 1):
				break
		if (models_counter > 1): # too many models
			print('>1 model')
			continue

		####
		print('OK')



		for residue in structure.get_residues():
			if (residue.get_id()[0].strip() != ''): # it's not a residue, may be HOH or H2SO4
				continue

			residue_name = residue.get_resname()
			residue_short_name = seq1(residue_name)

			if (residue_short_name == ''): # bad residue
				continue

			if residue_short_name in used: # we have already a sample of this amino acid
				continue

			print(residue)
			################ Saving the residue as PDBIO
			id = residue.get_id()[1]
			print('!!' + str(id))

			class Selector(Select):
				def __init__(self, id):
					self.id = id
					self.selected = False
				def accept_residue(self, residue):
					if residue.get_id()[1] == self.id and not self.selected:
						self.selected = True
						print(residue.get_resname() + '!!!!!!!!!!!!!!!!')
						return 1
					else:
						return 0
				def accept_atom(self, atom):
					if (not atom.is_disordered() or atom.get_altloc() == 'A'):
						atom.set_altloc(' ')
						return 1
					else:
						return 0

			writer.set_structure(structure)
			writer.save(os.path.join(result_folder, residue_name + '.ent'), Selector(id))
			################

			used.add(residue_short_name)
			if (len(used) == 20): # we have all the residues
				print('Done!')
				sys.exit(0)