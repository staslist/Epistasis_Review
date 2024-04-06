from ctypes import *
import random
from pysnptools.snpreader import Bed
import numpy as np
import pandas as pd
from scipy.stats import norm
import sys
import os

remma_path = os.path.split(os.path.realpath(__file__))[0]
fun_path = remma_path + '/fun'


def REMMA_epi_random(plink_prefix, output_prefix, p_val = 1.0/100000, num_test = 50000, num_random_pair = 1000000):
	
	snp_bed = Bed(plink_prefix, count_A1=False)
	num_ID = snp_bed.iid_count
	num_SNP = snp_bed.sid_count
	
	random_SNP_file = output_prefix + '.random_SNP'
	fout = open(random_SNP_file, 'w')
	fout.write('SNP1\tSNP2\n')
	SNP_vec = range(1, num_SNP+1)
	for i in range(num_random_pair):
		temp = random.sample(SNP_vec, 2)
		stri = str(temp[0]) + ' ' + str(temp[1]) + '\n'
		fout.write(stri)
	fout.close()
	
	C_path = fun_path + "/REMMA_epi_select_C"
	REMMA_epi_select_C = np.ctypeslib.load_library(C_path, ".")
	REMMA_epi_select_C.remma_epi_select_c.argtypes = [c_long, c_long, c_long, c_long, POINTER(c_char), POINTER(c_char), POINTER(c_char), c_double]
	REMMA_epi_select_C.remma_epi_select_c.restype = c_int
	
	x = REMMA_epi_select_C.remma_epi_select_c(num_ID, num_SNP, num_random_pair, num_test, plink_prefix, output_prefix, random_SNP_file, 0)
	
	random_SNP_res_file = random_SNP_file + '_res'
	data = pd.read_table(random_SNP_res_file, sep='\s+')
	
	data_mean = np.mean(data['epi_effect'])
	data_var = np.var(data['epi_effect'])
	quan_val = norm.isf(p_val, loc = 0, scale = np.sqrt(data_var))
	data_cor = np.corrcoef(np.abs(data['epi_effect']), data['chi_value'])
	return quan_val, data_mean, data_var, data_cor[0, 1]



if(sys.argv[1] == '--help'):
	print "Quick start: remma_epi_random --plink_prefix plink --out_prefix pre_file "
	print "--quantile  1e-5 --num_test 50000 --num_random_pair 1000000"
	sys.exit()

print "Read the parameters."
if(len(sys.argv) != 11):
	print "The program need 5 parameters! Please check!"
	sys.exit()

num_par = len(sys.argv)
label = 0
for i in range(1, 6):
	if(sys.argv[2*i-1] == '--plink_prefix'):
		plink_prefix = sys.argv[2*i]
		label = 1
if(label == 0):
	print "Please provide the parameter --plink_prefix"
	sys.exit()

label = 0
for i in range(1, 6):
	if(sys.argv[2*i-1] == '--out_prefix'):
		output_prefix = sys.argv[2*i]
		label = 1
if(label == 0):
	print "Please provide the parameter --out_prefix"
	sys.exit()


label = 0
for i in range(1, 6):
	if(sys.argv[2*i-1] == '--quantile'):
		p_val = float(sys.argv[2*i])
		label = 1
if(label == 0):
	print "Please provide the parameter --quantile"
	sys.exit()



label = 0
for i in range(1, 6):
	if(sys.argv[2*i-1] == '--num_test'):
		num_test = int(sys.argv[2*i])
		label = 1
if(label == 0):
	print "Please provide the parameter --num_test"
	sys.exit()


label = 0
for i in range(1, 6):
	if(sys.argv[2*i-1] == '--num_random_pair'):
		num_random_pair = int(sys.argv[2*i])
		label = 1
if(label == 0):
	print "Please provide the parameter --num_random_pair"
	sys.exit()



a, b, c, d = REMMA_epi_random(plink_prefix, output_prefix, p_val = p_val, num_test = num_test, num_random_pair = num_random_pair)

print "Information from the ", num_random_pair, " SNP pairs: "

print "The ", p_val, " quantile: ", a
print "The SNP-SNP epistatic sample mean and variance: ", b, c
print "The coefficient of the epistatic effects and chi-square value: ",d
