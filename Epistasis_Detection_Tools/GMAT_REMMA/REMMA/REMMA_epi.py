from ctypes import *
import os
import numpy as np
from pysnptools.snpreader import Bed
from pandas import DataFrame
import pandas as pd
from scipy.stats import chi2
import scipy.linalg
import sys


remma_path = os.path.split(os.path.realpath(__file__))[0]
fun_path = remma_path + '/fun'
sys.path.append(fun_path)
from start_end_snp import start_end_snp

def remma_epi(plink_prefix, output_prefix,
total_part = 1, which_part = 1, num_test_each = 50000, p_value = 0.0001):
	
	snp_bed = Bed(plink_prefix, count_A1=False)
	num_ID = snp_bed.iid_count
	num_SNP = snp_bed.sid_count
	
	start_SNP, end_SNP = start_end_snp(num_SNP, total_part, which_part)
	print start_SNP, end_SNP
	
	chi_value = chi2.isf(p_value, 1)
	
	C_path = fun_path + "/REMMA_epi_C"
	REMMA_epi_C = np.ctypeslib.load_library(C_path, ".")
	print REMMA_epi_C.remma_epi_c
	
	REMMA_epi_C.remma_epi_c.argtypes = [c_long, c_long, c_long, c_long, c_long, c_long, c_long, POINTER(c_char), POINTER(c_char), c_double]
	REMMA_epi_C.remma_epi_c.restype = c_int
	
	x = REMMA_epi_C.remma_epi_c(num_ID, num_SNP, start_SNP, end_SNP, total_part, which_part, num_test_each, plink_prefix, output_prefix, chi_value)
	



if(sys.argv[1] == '--help'):
	print "Quick start: remma_epi --plink_prefix plink --out_prefix pre_file "
	print "--parallel total i --num_test 50000 --p_threshold 1e-5"
	sys.exit()

print "Read the parameters."
if(len(sys.argv) != 12):
	print "The program need 5 parameters! Please check!"
	sys.exit()

num_par = len(sys.argv)

label = 0
for i in range(1, num_par):
	if(sys.argv[i] == '--plink_prefix'):
		plink_prefix = sys.argv[i+1]
		label = 1
if(label == 0):
	print "Please provide the parameter --plink_prefix"
	sys.exit()

label = 0
for i in range(1, num_par):
	if(sys.argv[i] == '--out_prefix'):
		output_prefix = sys.argv[i+1]
		label = 1
if(label == 0):
	print "Please provide the parameter --out_prefix"
	sys.exit()

label = 0
for i in range(1, num_par):
	if(sys.argv[i] == '--parallel'):
		total_part = long(sys.argv[i+1])
		which_part = long(sys.argv[i+2])
		label = 1
if(label == 0):
	print "Please provide the parameter --parallel"
	sys.exit()

label = 0
for i in range(1, num_par):
	if(sys.argv[i] == '--num_test'):
		num_test_each = long(sys.argv[i+1])
		label = 1
if(label == 0):
	print "Please provide the parameter --num_test"
	sys.exit()

label = 0
for i in range(1, num_par):
	if(sys.argv[i] == '--p_threshold'):
		p_threshold = float(sys.argv[i+1])
		label = 1
if(label == 0):
	print "Please provide the parameter --p_threshold"
	sys.exit()



remma_epi(plink_prefix, output_prefix,
total_part, which_part, num_test_each = num_test_each, p_value = p_threshold)
