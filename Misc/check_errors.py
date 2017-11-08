import os


"""

Checks if there are no wrong kmers

Boris Starkov(c)

"""

data_folder = '../Data'

kmers_path = os.path.join(data_folder, 'fixed_size_chains')

files_processed = 0
with open(kmers_path, 'r') as file:
	while True:
		seq = file.readline().strip()
		if seq == '':
			break

		#
		files_processed += 1
		if (files_processed % 100000 == 0):
			print('processed %d files' % files_processed)
		#
		angles_string = file.readline().strip()
		angles = angles_string.split(' ')

		if len(angles) != 3 * (len(seq) - 1):
			print(seq + ' is bad!')
			print(str(len(seq)) + " " + str(len(angles)))
			print(seq)
			print('!'.join(angles))
			sys.exit(1)

		for angle in angles:
			try:
				float(angle)
			except ValueError:
				print(seq + ' is bad! and has non-float value')
				sys.exit(1)

print('Done! Alright.')