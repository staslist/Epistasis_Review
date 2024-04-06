import numpy as np
import pandas as pd
from scipy import linalg
import common
from scipy.optimize import minimize

def opt_EMMA_sigma(y, kin, X = None, init = [1.0, 1.0], lim = [10e-10, 10e10]):
	
	if X is None:
		print 'Fixed design matrix is not provided, intercept is included!'
		X = np.array([1.0]*len(y)).reshape(len(y),1)
	
	D, U = linalg.eigh(kin) #eigen decomposition kinship
	D.shape = 1, len(y)
	
	y = np.array(y)
	y.shape = len(y),1
	Xrot = np.dot(U.T, X) #rotate 
	yrot = np.dot(U.T, y)
	
	###-2L
	def LLF(var):
		
		addVar,resVar = var
		
		#log|V|
		V_diag = D*addVar + resVar
		LL = np.sum(np.log(V_diag))
		#print LL
		
		#log|X'VX|
		Vinv = 1.0/V_diag
		XVX = common.Dtri_matT(Xrot.T, Vinv)
		LL += np.linalg.slogdet(XVX)[1]
		#print LL
		
		#y'Py
		LL += common.Dtri_matT(yrot.T, Vinv)[0,0] #y'Vy
		#print LL
		yVX = common.Dtri_mat(yrot.T, Vinv, Xrot)
		XVXinv = np.linalg.inv(XVX)
		LL -= common.tri_matT(yVX, XVXinv)[0,0]
		#print LL
		
		return LL
	lim = [lim]*2
	result = minimize(LLF, init, method = 'L-BFGS-B', bounds = lim)
	
	return result 
