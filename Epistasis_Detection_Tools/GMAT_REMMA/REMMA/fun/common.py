import numpy as np
from collections import defaultdict
def is_int(num):
	try:
		int(num)
		return True
	except ValueError:
		return False

def is_float(num):
	try:
		float(num)
		return True
	except ValueError:
		return False

###declare the 2D dict
def dct_3D():
	return defaultdict(dct_31D)

def dct_31D():
	return defaultdict(dct_32D)

def dct_32D():
	return defaultdict()


###declare the 2D dict
def dct_2D():
	return defaultdict(dct_21D)

def dct_21D():
	return defaultdict()


###declare the 1D dict
def dct_1D():
	return defaultdict()

###triple matrix multiplication a x b x a'

def tri_matT(a, b):
	res = np.dot(a,b)
	res = np.dot(res,a.T)
	return res

###triple matrix multiplication a x b x c

def tri_mat(a, b, c):
	res = np.dot(a, b)
	res = np.dot(res, c)
	return res

###triple matrix multiplication a x D x a'  D is diagonal matrix in row vector form

def Dtri_matT(a, b):
	res = np.multiply(a,b)
	res = np.dot(res,a.T)
	return res

###triple matrix multiplication a x D x c  D is diagonal matrix in row vector form

def Dtri_mat(a, b, c):
	res = np.multiply(a,b)
	res = np.dot(res, c)
	return res

