import os
import sys
import numpy as np
from pysnptools.snpreader import Bed
from pandas import DataFrame
import pandas as pd
from scipy.stats import chi2
import gc

remma_path = os.path.split(os.path.realpath(__file__))[0]
fun_path = remma_path + '/fun'
sys.path.append(fun_path)
from opt_EMMA_sigma import opt_EMMA_sigma
from addNormMat import addNormMat

def remma_add(plink_prefix, pheno_file, out_file):
	
	print "Read the phenotypic file"
	pheno = pd.read_table(pheno_file, header = None)
	snp_bed = Bed(plink_prefix, count_A1 = False)
	num_ID = snp_bed.iid_count
	
	if(pheno.shape[0] != num_ID):
		print "Individual number in phenotypic file is not equal to that in plink file! Please check!"
		sys.exit()
	
	if(pheno.shape[1] < 4):
		print "Please provide at least 4 columns (family ID, individual ID, intercept and phenotypic values) in the phenotypic file!"
		sys.exit()
	
	
	X = np.array(pheno.ix[:,2:(pheno.shape[1]-2)]).astype(np.float)
	X.shape = num_ID, pheno.shape[1] - 3
	
	y = np.array(pheno.ix[:,pheno.shape[1]-1]).astype(np.float)
	y.shape = num_ID, 1
	
	
	print "Read the kinship file:"
	add_kin_file = plink_prefix + '.addGmat.grm0'
	add_kin = np.loadtxt(add_kin_file)
	
	print "Estimate the varinaces"
	val = opt_EMMA_sigma(y, add_kin, X)
	if(not val['success']):
		print 'Not Converge for Variance!!!'
		sys.exit()
	
	add_var, res_var = val['x']
	print "The addtive and residual variances:", add_var, res_var
	
	
	print "Build the coef of MME"
	coef_top = np.dot(X.T, X)
	coef_top = np.concatenate((coef_top, X.T), axis=1)
	
	add_kin_inv_file = plink_prefix + '.addGmat.giv0'
	add_kin_inv = np.loadtxt(add_kin_inv_file)
	coef_bottom = np.eye(num_ID)
	coef_bottom = np.add(np.multiply(add_kin_inv, res_var/add_var), coef_bottom)
	coef_bottom = np.concatenate((X, coef_bottom), axis=1)
	
	coef = np.concatenate((coef_top, coef_bottom), axis=0)
	del coef_top,coef_bottom
	gc.collect()
	coef = np.linalg.inv(coef)
	
	
	print "Build the right hand of MME"
	right_hand_up = np.dot(X.T, y)
	right_hand = np.concatenate((right_hand_up, y), axis=0)
	
	print "Cal the additve effect and PEV"
	effect = np.dot(coef, right_hand)
	add_effect = effect[(-num_ID):,:]
	kin_inv_by_effect = np.dot(add_kin_inv, add_effect)
	add_effect_var = np.subtract(np.multiply(add_kin, add_var), np.multiply(coef[(-num_ID):, (-num_ID):], res_var))
	snp_effect_var_mid = np.dot(np.dot(add_kin_inv, add_effect_var), add_kin_inv)
	del coef, add_kin, add_kin_inv
	gc.collect()
	
	
	print "Read the SNP matrix"
	snp_norm = addNormMat(plink_prefix)
	
	print "Association"
	add_snp_effect = np.dot(snp_norm, kin_inv_by_effect)
	add_snp_effect = add_snp_effect[:,0]
	
	snp_norm = np.dot(snp_norm, snp_effect_var_mid) * snp_norm
	add_snp_var = np.sum(snp_norm, axis=1)
	add_snp_chi = add_snp_effect*add_snp_effect/add_snp_var
	add_snp_p = chi2.sf(add_snp_chi, 1)
	
	result = {'SNPID': snp_bed.sid,
		'chro': snp_bed.pos[:,0],
		'chrPos': snp_bed.pos[:,2],
		'effect': add_snp_effect,
		'chi_val': add_snp_chi,
		'Pvalue': add_snp_p}
	
	resultFrame = DataFrame(result, columns=['SNPID', 'chro', 'chrPos', 'effect', 'chi_val', 'Pvalue'])
	
	resultFrame.to_csv(out_file, index = False, header = True, sep = '\t')



if(sys.argv[1] == '--help'):
	print "Quick start: remma_add --prefix plink --pheno phe --out out_file"
	sys.exit()

print "Read the parameters."
if(len(sys.argv) != 7):
	print "The program need 3 parameters! Please check!"
	sys.exit()

label = 0
for i in range(1, 4):
	if(sys.argv[2*i-1] == '--prefix'):
		plink_prefix = sys.argv[2*i]
		label = 1
if(label == 0):
	print "Please provide the parameter --prefix"
	sys.exit()

label = 0
for i in range(1, 4):
	if(sys.argv[2*i-1] == '--pheno'):
		pheno_file = sys.argv[2*i]
		label = 1
if(label == 0):
	print "Please provide the parameter --pheno"
	sys.exit()

label = 0
for i in range(1, 4):
	if(sys.argv[2*i-1] == '--out'):
		out_file = sys.argv[2*i]
		label = 1
if(label == 0):
	print "Please provide the parameter --out"
	sys.exit()


remma_add(plink_prefix, pheno_file, out_file)
