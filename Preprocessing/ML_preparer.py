import os
import time

import math
import random

"""

Prepares the k-mers to machine learning:
*) Divides them into 3 file : train, test and closed_test
*) Divides every angle by PI so we could predict [-1;1] numbers intead of [-PI;PI]

Boris Starkov(c)

"""
ninemers_file = 'fixed_size_chains'

objects_to_process = -1 # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ very important param. How much objects to take. if = -1 all the objects will be taken.

processed_objects = 0




data_folder = '../Data'

kmers_path = os.path.join(data_folder, ninemers_file)

#
datasets_folder = os.path.join(data_folder, 'Datasets')
if not os.path.exists(datasets_folder):
	os.makedirs(datasets_folder)
#


train_path = os.path.join(datasets_folder, 'train')
test_path = os.path.join(datasets_folder, 'test')
closed_test_path = os.path.join(datasets_folder, 'closed_test')

train_printer = open(train_path, 'w')
test_printer = open(test_path, 'w')
closed_test_printer = open(closed_test_path, 'w')

train_probability = 70 # percent of tests that go to train file
test_probability = 25 # percent of tests that go to test file
# others go to closed test

part = 1 # percantage of kmers to pass

def work(s): #makes float from strign and divides it by pi
	res = float(s) / math.pi
	res = min(res, 1)
	res = max(res, -1)
	return res

start_time = time.time()

#
print('Started dividing kmers into samples')
#

with open(kmers_path) as file:
	while True:
		sequence_line = file.readline().rstrip()
		if (sequence_line == ''):
			break
		
		angles_line = file.readline().rstrip()
		angles = map(work, angles_line.split(' '))

		random_integer = random.randint(0, 99)

		if random_integer >= part:
			continue
		
		random_integer = random.randint(0, 99)

		if random_integer < train_probability:
			train_printer.write(sequence_line + '\n')
			train_printer.write(' '.join(map(str, angles)) + '\n')

		elif random_integer < (train_probability + test_probability):
			test_printer.write(sequence_line + '\n')
			test_printer.write(' '.join(map(str, angles)) + '\n')

		else:
			closed_test_printer.write(sequence_line + '\n')
			closed_test_printer.write(' '.join(map(str, angles)) + '\n')

		processed_objects  += 1
		print(processed_objects)
		if processed_objects == objects_to_process:
			break

print('Done in ' + str(time.time() - start_time) + ' seconds.')

train_printer.close()
test_printer.close()
closed_test_printer.close()