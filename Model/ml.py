import numpy as np
import keras

from keras import optimizers

from keras.models import Sequential
from keras.layers import Dense, Activation, Flatten, Dropout

import os

import math
import random

k = 9 # number of residues
in_k = 22 * k # size of input matrix
out_k = 3 * (k - 1) # number of angles

X_train = np.array([])
y_train = np.array([])

X_test = np.array([])
y_test = np.array([])


train_filename = '../Data/Datasets/train'
test_filename = '../Data/Datasets/test'

model_filename = 'ML_model.h5'


dict = {'a':0, 'r':1, 'n':2,
'd':3, 'c':4, 'e':5, 'q':6,
'g':7, 'h':8, 'i':9, 'l':10,
'k':11, 'm':12, 'f':13, 'p':14,
's':15, 't':16, 'w':17,
'y':18, 'v':19, 'u':20, 'o':21}

def char_to_vector(c):
	res = [0] * 22
	one_pos = dict[c]
	res[one_pos] = 1
	return res

def seq_to_matrix(s):
	res = []
	for c in s.lower():
		res.append(char_to_vector(c))
	return res


def fill_data(filename): # fills with data from file
	
	#
	lines_counter = 0
	#

	X_matrixes = []
	y_vectors = []

	file = open(filename, 'r')

	while True:
		sequence = file.readline().lower().strip()
		if (sequence == ''):
			break
		angles_string = file.readline().strip()

		#
		lines_counter += 1
		if lines_counter % 1000 == 0:
			print('%d lines processed' % lines_counter)
		#
		#print(sequence + '!')
		angles = map(float, angles_string.split(' '))

		X_matrixes.append(seq_to_matrix(sequence))
		y_vectors.append(angles)

	""" WARNING: old version
	lines = file.readlines() # CAUTION! it loads the whole file so we suppose file is not so large
	#random.shuffle(lines) # model.fit does it

	for line in lines:
		line = line.replace('\n', '')
		
		lines_counter += 1
		if (lines_counter % 1000 == 0):
			print('%d lines processed' % lines_counter)
		#

		data = line.split('\t')
		seq = data[0]
		angles = [float(x) for x in data[1:]]

		X_matrixes.append(seq_to_matrix(seq))
		y_vectors.append(angles)
	"""

	if (filename == train_filename):
		global X_train, y_train
		X_train = np.array(X_matrixes)
		y_train = np.array(y_vectors)
		
		#
		print(X_train)
		print('///')
		print(y_train)
		#
		
	else:
		global X_test, y_test
		X_test = np.array(X_matrixes)
		y_test = np.array(y_vectors)

		#
		print(X_test)
		print('///')
		print(y_test)
		#


activation = 'tanh'
epochs = 50



fill_data(test_filename)

#ML starts here

#
def load():
	if not os.path.exists(model_filename):
		return False
	res = raw_input('Do you want to load prev model? (y/n)').lower()
	return res == 'y'

if load():
	
	#
	print('the model was loaded')
	#

	model = keras.models.load_model(model_filename)
else:

	#
	print('the model trains again')
	#

	fill_data(train_filename)
	
	model = Sequential()
	#r = (float(in_k) / out_k) ** (float(1) / 3)
	neurons_1 = 96 # out * r * r = 96
	neurons_2 = 48 # out * r = 48

	model.add(Flatten(input_shape = (k, 22)))
	model.add(Dense(activation = activation, use_bias = True, units=neurons_1 )) # 1st hidden layer
	model.add(Dropout(0.25))
	model.add(Dense(activation = activation, use_bias = True, units = neurons_2 )) # 2nd hidden layer
	model.add(Dropout(0.25))
	model.add(Dense(activation = activation, use_bias = True, units = out_k)) # output layer


	model.compile(loss='mean_squared_error', metrics=['mse', 'mae'], optimizer=optimizers.SGD(lr=0.1, decay=1e-6, momentum=0.95, nesterov=True))

	model.fit(x = X_train, y = y_train, epochs = epochs, batch_size = 128)

model.save(model_filename)

print('Test set: ')
print(model.evaluate(x = X_test, y = y_test, batch_size = 128))