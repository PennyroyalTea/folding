import os
import time
import sys
"""

Takes proteins of length 9 and writes them as cutter.py does.

Boris Starkov(c)

"""

k = 9 # size of new sequences

data_folder = '../Data'

processed_files_folder = os.path.join(data_folder, 'processed_pdb')

result_path = os.path.join(data_folder, 'fixed_size_chains_original')

printer = open(result_path, 'w')

start_time = time.time()

##
print('started cutting')
##

for file_name in os.listdir(processed_files_folder):
	##
	print('now working on ' + file_name)
	##
	file_path = os.path.join(processed_files_folder, file_name)

	with open(file_path) as file:
		content = file.read()

	lines = content.split('\n')

	for i in range(0, len(lines) - 1, 2):
		sequence = lines[i].strip()
		angles = lines[i + 1].strip().split(' ')

		if len(sequence) != k:
			continue
		
		printer.write(sequence + '\n')
		printer.write(' '.join(angles) + '\n')

##
print('Done in ' + str(time.time() - start_time) + ' seconds!')		
##
printer.close()