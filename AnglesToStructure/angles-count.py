import os
import sys
import gzip

from Bio.PDB import *
from Bio import SeqIO
from Bio.SeqUtils import seq1

import warnings

from collections import defaultdict

"""

counts angles NCaC in every aminoacid

Boris Starkov(c)

"""

RESIDUES_TO_PROCESS = 100000

residues_processed = 0

angles = defaultdict(list)

good_residues = {'ala', 'arg', 'asn', 'asp', 'cys', 'glu', 'gln', 'gly', 'his', 'ile', 'leu', 'lys', 'met', 'phe', 'pro', 'ser', 'thr', 'trp', 'tyr', 'val'}




warnings.filterwarnings('ignore')

data_folder = '../Data'
PDB_files_directory = os.path.join(data_folder, 'PDB_files')

result_folder = os.path.join(data_folder, 'Residue_samples')
if not os.path.exists(result_folder):
	os.makedirs(result_folder)


parser = PDBParser()

BREAK = False

#go through all the files
for root, dirs, files in os.walk(PDB_files_directory):
	
	if (BREAK): 
		break

	for file_name in files:

		if (BREAK):
			break

		file_path = os.path.join(root, file_name) # path to current file

		with gzip.open(file_path, 'rb') as file:
			structure = parser.get_structure('', file)

		#### checks if a structure fits


		## bad resolution
		if (structure.header['resolution'] == None or structure.header['resolution'] > 3): #bad resolution
			print(file_name + ' skipped: bad resolution ' + str(structure.header['resolution']))
			continue 

		## more than 1 model
		models_counter = 0
		for model in structure.get_models():
			models_counter += 1
			if (models_counter > 1):
				break
		if (models_counter > 1): # too many models
			print(file_name + ' skipped: more than one model')
			continue

		####
		print(file_name + ' is OK')



		for residue in structure.get_residues():
			if (residue.get_id()[0].strip() != ''): # it's not a residue, may be HOH or H2SO4
				continue

			residue_name = residue.get_resname()
			residue_short_name = seq1(residue_name)

			if not residue_name.lower() in good_residues: # bad residue
				continue

			
			atoms = []

			has_disordered = False

			for atom in residue.get_atoms():
				if (atom.is_disordered() or atom.get_altloc() == 'A'):
					has_disordered = True
					break
				atoms.append(atom.get_vector())
				if (len(atoms) == 3):
					break

			if (len(atoms) < 3): # less than 3 atoms in a residue
				continue

			if (has_disordered):
				continue

			cur_angle = calc_angle(atoms[0], atoms[1], atoms[2])

			angles[residue_name].append(cur_angle)

			residues_processed += 1
			print(residues_processed)
		if (residues_processed >= RESIDUES_TO_PROCESS):
			BREAK = True
			break

mean_overall = 0

for name, cur_angles in angles.items():
	num = len(cur_angles)
	mean = 1.0 * sum(cur_angles) / len(cur_angles)
	mean_overall += mean
	mst = 0.0
	for cur_angle in cur_angles:
		mst += (mean - cur_angle) ** 2
	mst /= num
	mst = mst ** (0.5)
	print(name + ' ' + str(mean - mst) + ' - ' + str(mean + mst))

mean_overall /= 20
print('mean overall is ' + str(mean_overall))