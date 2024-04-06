import pandas as pd
from pysnptools.snpreader import Bed
import sys
import numpy as np
import gc
import scipy
from scipy.stats import chi2
from pandas import DataFrame
import os

remma_path = os.path.split(os.path.realpath(__file__))[0]
fun_path = remma_path + '/fun'
sys.path.append(fun_path)
import opt_multikin
from addNormMat import addNormMat



def remma_pre(plink_prefix, pheno_file, output_prefix, add_test = 0):
	
	print "Read the phenotypic file"
	pheno = pd.read_table(pheno_file, header = None)
	snp_bed = Bed(plink_prefix, count_A1 = False)
	
	if(pheno.shape[0] != snp_bed.iid_count):
		print "Individual number in phenotypic file is not equal to that in plink file! Please check!"
		sys.exit()
	
	if(pheno.shape[1] < 4):
		print "Please provide at least 4 columns (family ID, individual ID, intercept and phenotypic values) in the phenotypic file!"
		sys.exit()
	
	X = np.array(pheno.ix[:,2:(pheno.shape[1]-2)]).astype(np.float)
	X.shape = snp_bed.iid_count, pheno.shape[1] - 3
	y = np.array(pheno.ix[:,pheno.shape[1]-1]).astype(np.float)
	y.shape = snp_bed.iid_count, 1
	
	print "Read the kinship file"
	add_kin_file = plink_prefix + '.addGmat.grm0'
	epi_kin_file = plink_prefix + '.epiGmatAA.grm0'
	
	add_kin = np.loadtxt(add_kin_file)
	epi_kin = np.loadtxt(epi_kin_file)
	
	kin = [add_kin,epi_kin]
	
	
	print "Estimate the variances"
	val = opt_multikin.opt_multikin(y, kin, X)
	if(not val['success']):
		print "Not Converge for Variance!!!";
		sys.exit()
	
	add_var, epi_var, res_var = val['x'][0], val['x'][1], val['x'][2]
	print "The additive, epistatic and residual variances are: ", val['x'][0], val['x'][1], val['x'][2]
	
	
	print 'Build the coef matrix for MME'
	
	#up
	coef_up = np.dot(X.T, X)
	coef_up = np.concatenate((coef_up, X.T, X.T), axis=1)
	
	#middle
	add_kin_inv_file = plink_prefix + '.addGmat.giv0'
	add_kin_inv = np.loadtxt(add_kin_inv_file)
	diag_one = np.eye(snp_bed.iid_count)
	coef_mid = np.add(add_kin_inv*res_var/add_var, diag_one)
	coef_mid = np.concatenate((X, coef_mid, diag_one), axis=1)
	
	#bottom
	epi_kin_inv_file = plink_prefix + '.epiGmatAA.giv0'
	epi_kin_inv = np.loadtxt(epi_kin_inv_file)
	coef_bottom = np.add(epi_kin_inv*res_var/epi_var, diag_one)
	coef_bottom = np.concatenate((X, diag_one, coef_bottom), axis=1)
	del diag_one
	gc.collect()
	
	coef = np.concatenate((coef_up, coef_mid, coef_bottom), axis=0)
	del coef_up, coef_mid, coef_bottom
	gc.collect()
	
	#all
	coef = np.linalg.inv(coef)
	
	
	print "Build the right hand of MME"
	right_hand = np.dot(X.T, y)
	right_hand = np.concatenate((right_hand, y, y), axis=0)
	
	#effect
	num_ID = snp_bed.iid_count
	effect = np.dot(coef, right_hand)
	
	#test the additive effect
	if(addTest==1):
		print "Test the additive effect"
		add_effect = effect[(-2*num_ID):(-num_ID), 0]
		kin_inv_by_effect = np.dot(add_kin_inv, add_effect)
		addMarkerMat = addNormMat(plink_prefix)
		snp_effect = np.dot(addMarkerMat, kin_inv_by_effect)
		add_effect_var = np.subtract(add_kin*add_var, coef[(-2*num_ID):(-num_ID), (-2*num_ID):(-num_ID)]*res_var)
		add_effect_var = np.dot(np.dot(add_kin_inv, add_effect_var), add_kin_inv)
		add_snp_var = np.multiply(np.dot(addMarkerMat, add_effect_var), addMarkerMat)
		del addMarkerMat, add_effect_var
		gc.collect()
		add_snp_var = np.sum(add_snp_var, axis=1)
		add_snp_chi = np.divide(np.multiply(snp_effect, snp_effect), add_snp_var)
		
		add_snp_p = chi2.sf(add_snp_chi, 1)
		
		result = {'SNPID': snp_bed.sid,
			'chro': snp_bed.pos[:,0],
			'chrPos': snp_bed.pos[:,2],
			'effect': snp_effect,
			'chi_val': add_snp_chi,
			'Pvalue': add_snp_p}
		
		resultFrame = DataFrame(result, columns=['SNPID', 'chro', 'chrPos', 'effect', 'chi_val', 'Pvalue'])
		add_output_file = out_file + '.addres'
		resultFrame.to_csv(add_output_file, index = False, header = True, sep = ' ')
	
	del add_kin, add_kin_inv
	gc.collect()
	
	print "Prepare the temp file for the next step"
	
	#effect file
	epi_effect = effect[(-num_ID):,:]
	kin_inv_by_effect = np.dot(epi_kin_inv, epi_effect)
	kin_inv_by_effect_file = out_file + '.eff'
	np.savetxt(kin_inv_by_effect_file, kin_inv_by_effect)
	
	#error file
	epi_effect_var = np.subtract(epi_kin*epi_var, coef[(-num_ID):, (-num_ID):]*res_var)
	del coef,epi_kin
	gc.collect()
	
	epi_effect_var = np.dot(np.dot(epi_kin_inv, epi_effect_var), epi_kin_inv)
	del epi_kin_inv
	gc.collect()
	
	D, U = scipy.linalg.eigh(epi_effect_var)
	
	eigenvec_file = out_file + '.eigen_vec'
	eigenval_file = out_file + '.eigen_val'
	np.savetxt(eigenvec_file, U)
	np.savetxt(eigenval_file, D)

if(sys.argv[1] == '--help'):
	print "Quick start: remma --plink_prefix plink --pheno phe --out_prefix out_file --add_test 0/1"
	sys.exit()
print "Read the parameters"
if(len(sys.argv) != 9):
	print "The program need 4 parameters! Please check!"
	sys.exit()

label = 0
for i in range(1, 5):
	if(sys.argv[2*i-1] == '--plink_prefix'):
		plink_prefix = sys.argv[2*i]
		label = 1
if(label == 0):
	print "Please provide the parameter --plink_prefix"
	sys.exit()

label = 0
for i in range(1, 5):
	if(sys.argv[2*i-1] == '--pheno'):
		pheno_file = sys.argv[2*i]
		label = 1
if(label == 0):
	print "Please provide the parameter --pheno"
	sys.exit()

label = 0
for i in range(1, 5):
	if(sys.argv[2*i-1] == '--out_prefix'):
		out_file = sys.argv[2*i]
		label = 1
if(label == 0):
	print "Please provide the parameter --out_prefix"
	sys.exit()

label = 0
for i in range(1, 5):
	if(sys.argv[2*i-1] == '--add_test'):
		addTest = int(sys.argv[2*i])
		label = 1
if(label == 0):
	print "Please provide the parameter --add_test"
	sys.exit()


remma_pre(plink_prefix, pheno_file, out_file, add_test = addTest)
