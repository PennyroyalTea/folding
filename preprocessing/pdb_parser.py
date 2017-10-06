import random
import warnings
import time

import os
import gzip

import multiprocessing
from multiprocessing import Process

from Bio.PDB import *
from Bio import SeqIO
from Bio.SeqUtils import seq1

"""

Unzippes, filter and parses .ent files to one big .txt with residues of length 9.

Boris Starkov(c)

"""

warnings.filterwarnings('ignore')

data_folder = '../Data'

PDB_files_folder = os.path.join(data_folder, 'PDB_files')

PDB_processed_files_folder = os.path.join(data_folder, 'processed_pdb')

if not os.path.exists(PDB_processed_files_folder):
	os.makedirs(PDB_processed_files_folder)

result_file = os.path.join(PDB_processed_files_folder, 'data')

log_folder = os.path.join(data_folder, 'logs')
log_file = os.path.join(log_folder, 'parser')


n_CPU = multiprocessing.cpu_count() # number of CPU on your PC. May be found out using lscpu command.
print(str(n_CPU) + ' threads in work.')

parser = PDBParser()


random.seed(566)

def work(folders_to_process, id):
	counter = 0

	global printer
	printer = open(result_file + '_' + str(id), 'w')

	global logger
	logger = open(log_file + '_' + str(id) + '.log', 'w')

	for folder_path in folders_to_process:
		for file_name in os.listdir(folder_path):
			counter += 1
			file_path = os.path.join(folder_path, file_name) # full path to the file
			
			with gzip.open(file_path, 'rb') as finput:
				structure = parser.get_structure('', finput)


			#
			logger.write('#' + str(counter) + ' ' + file_name[3:7] + '\n')
			#

			handle_PDB(structure, file_name[3:7])

	printer.close()
	logger.close()

def check(structure): # checks if a structure is suitable for work
	"""
	0: structure is ok

	reasons to filter out the structure:
	1: bad resolution
	2: more than one model
	3: it's RNA/DNA or whatever
	4: it has disordered atoms
	"""


	#### checking the resolution
	resolution = structure.header['resolution']
	if (resolution == None or resolution > 3):
		logger.write('1 : bad resolution(' + str(resolution) + ') -> skipped\n')
		return False
	####

	#### checking if there is only one model
	model_counter = 0
	for model in structure.get_models():
		model_counter += 1
		if model_counter > 1:
			break
	if model_counter > 1:
		logger.write('2 : more than one model -> skipped\n')
		return False
	####

	#### check if it is a protein
	for residue in structure.get_residues():
		residue_name = residue.get_resname().lower().strip()
		if residue.get_id()[0].strip() == '' and seq1(residue_name) == '': # it's not a hetatm(water/acid) but also not a known aminoacid... what is it? probably a DNA!
			logger.write('3 : has unknown atoms(probably a DNA or RNA or whatever) -> skipped\n')
			return False
	####

	#### check if it has disordered atoms
	for atom in structure.get_atoms():
		if (atom.is_disordered()):
			logger.write('4 : one of atoms is disordered -> skipped\n')
			return False
	####
	
	logger.write('0 : structure is OK\n')
	return True


def handle_PDB(structure, name):
	if not check(structure): # structure is brokes some way
		return

	for chain in structure.get_chains():
		atoms = [] # list of main atoms: N Ca C
		sequence = [] # sequance of a protein

		for residue in chain.get_residues():

			residue_name = residue.get_resname().lower().strip() # name of residue in 3-letters style, e.g. gly

			if (residue.get_id()[0].strip() == ''): # if it's not a hetatm 
				sequence.append(seq1(residue_name)) # adds name of residue in 1-letter style

				atoms_added = 0
				for atom in residue.get_atoms():
					atoms.append(atom)
					atoms_added += 1
					if atoms_added == 3: # we are only interested in first 3 atoms : N Ca C
						break
		"""

		calculating angles from atom positions

		"""

		angles = []

		for i in range(0, len(atoms) - 5, 3):
			vec0 = atoms[i].get_vector()
			vec1 = atoms[i + 1].get_vector()
			vec2 = atoms[i + 2].get_vector()
			vec3 = atoms[i + 3].get_vector()
			vec4 = atoms[i + 4].get_vector()
			vec5 = atoms[i + 5].get_vector()



			angleA = float(calc_dihedral(vec0, vec1, vec2, vec3))
			angleB = float(calc_dihedral(vec1, vec2, vec3, vec4))
			angleC = float(calc_dihedral(vec2, vec3, vec4, vec5))

			angles.append(str(angleA))
			angles.append(str(angleB))
			angles.append(str(angleC))

		printer.write(''.join(sequence) + '\n')
		printer.write(' '.join(angles) + '\n')





if __name__ == '__main__': #it's a main thread
	begin_time = time.time()


	all_PDB_subfolders = os.listdir(PDB_files_folder)

	worklist = [] # worklist[i] - names of folders to be processed by i-th thread

	for i in range(n_CPU):
		worklist.append([])

	for subfolder_name in all_PDB_subfolders:
		full_path = os.path.join(PDB_files_folder, subfolder_name) # full path to the subfolder
		worklist[random.randint(0, n_CPU - 1)].append(full_path)

	processes = [0] * n_CPU
	for i in range(n_CPU):
		processes[i] = Process(target=work, args=(worklist[i], i)) # initialization of every process
		processes[i].start()
	for i in range(n_CPU):
		processes[i].join() # waiting for process to end
	
	print('Parsing done! time: ' + str(time.time() - begin_time) + ' seconds.')