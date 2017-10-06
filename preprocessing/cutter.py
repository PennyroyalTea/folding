import os
import time

"""

Cuts all the chains into chains of length k (usually, 9) and makes ONE file with them.

Boris Starkov(c)

"""

k = 9 # size of new sequences

data_folder = '../Data'

processed_files_folder = os.path.join(data_folder, 'processed_pdb')

result_path = os.path.join(data_folder, 'fixed_size_chains')

printer = open(result_path, 'w')

start_time = time.time()

for file_name in os.listdir(processed_files_folder):
	##
	print('now working on ' + file_name)
	##
	file_path = os.path.join(processed_files_folder, file_name)

	with open(file_path) as file:
		content = file.read()

	lines = content.split('\n')

	for i in range(0, len(lines) - 1, 2):
		sequence = lines[i]
		angles = lines[i + 1].split(' ')

		for j in range(len(sequence) - k + 1):
			cur_sequence = sequence[j : (j + k)]
			cur_angles = angles[(3 * j) : (3 * (j + k - 1))]

			printer.write(cur_sequence + '\n')
			printer.write(' '.join(cur_angles) + '\n')

print('Done in ' + str(time.time() - start_time) + ' seconds!')		

printer.close()