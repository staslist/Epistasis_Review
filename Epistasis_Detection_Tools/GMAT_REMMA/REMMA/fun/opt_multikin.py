import numpy as np
import pandas as pd
from scipy import linalg
import common
from scipy.optimize import minimize

def opt_multikin(y, kin, X = None, init = None, lim = [10e-10, 10e10]):
	
	N = len(y)
	num_var = (len(kin)+1)
	if X is None:
		print 'Fixed design matrix is not provided, intercept is included!'
		X = np.array([1.0]*N).reshape(N, 1)
	
	y = np.array(y)
	y.shape = N,1
	
	if init is None:
		print 'Initial values for variances are not provied, 1.0s are used!'
		init = [1.0]*num_var
	
	###-2L
	def LLF(var):
		
		#log|V|
		V = np.diag([var[-1]]*N)
		for i in range(len(kin)):
			V += np.multiply(kin[i], var[i])
		LL = np.linalg.slogdet(V)[1]
		#print LL
		
		#log|X'VX|
		V = np.linalg.inv(V)
		XVX = common.tri_matT(X.T, V)
		LL += np.linalg.slogdet(XVX)[1]
		#print LL
		
		#y'Py
		LL += common.tri_matT(y.T, V)[0,0] #y'Vy
		#print LL
		yVX = common.tri_mat(y.T, V, X)
		XVX = np.linalg.inv(XVX)
		LL -= common.tri_matT(yVX, XVX)[0,0]
		#print LL
		
		return LL
	lim = [lim]*num_var
	result = minimize(LLF, init, method = 'L-BFGS-B', bounds = lim)
	#result = minimize(LLF, init, method = 'BFGS')
	return result

