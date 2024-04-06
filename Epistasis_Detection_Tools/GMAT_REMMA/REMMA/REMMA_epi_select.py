from ctypes import *
import sys
import numpy as np
from pysnptools.snpreader import Bed
from scipy.stats import chi2
import os


remma_path = os.path.split(os.path.realpath(__file__))[0]
fun_path = remma_path + '/fun'

def remma_epi_select(plink_prefix, output_prefix, SNP_select_file,
num_test = 50000, p_value = 0.0001):
	
	snp_bed = Bed(plink_prefix, count_A1=False)
	num_ID = snp_bed.iid_count
	num_SNP = snp_bed.sid_count
	
	fin = open(SNP_select_file, "r")
	line = fin.readline()
	num_SNP_select = -1;
	while line:
		num_SNP_select += 1;
		line = fin.readline()
	
	fin.close()
	
	
	
	chi_value = chi2.isf(p_value, 1)
	
	C_path = fun_path + "/REMMA_epi_select_C"
	REMMA_epi_select_C = np.ctypeslib.load_library(C_path, ".")
	print REMMA_epi_select_C.remma_epi_select_c
	
	REMMA_epi_select_C.remma_epi_select_c.argtypes = [c_long, c_long, c_long, c_long,
	POINTER(c_char), POINTER(c_char), POINTER(c_char), c_double]
	REMMA_epi_select_C.remma_epi_select_c.restype = c_int
	
	
	x = REMMA_epi_select_C.remma_epi_select_c(num_ID, num_SNP, num_SNP_select, num_test,
	plink_prefix, output_prefix, SNP_select_file, chi_value)
	
	out_file = SNP_select_file + '_res'
	out_file2 = out_file + ".addP"
	
	fin = open(out_file, 'r')
	fout = open(out_file2, 'w')
	line = fin.readline()
	stri = line.strip()
	stri = stri + " p_value\n";
	fout.write(stri)
	line = fin.readline()
	while line:
		arr = line.split()
		p_val = chi2.sf(float(arr[-1]), 1)
		stri = line.strip()
		stri = stri + " " + str(p_val) + "\n"
		fout.write(stri)
		line = fin.readline()
	fin.close()
	fout.close()



if(sys.argv[1] == '--help'):
	print "Quick start: remma_epi_select --plink_prefix plink --out_prefix pre_file "
	print "--SNP_select_file file --num_test 50000 --p_val 1"
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
	if(sys.argv[2*i-1] == '--SNP_select_file'):
		SNP_select_file = sys.argv[2*i]
		label = 1
if(label == 0):
	print "Please provide the parameter --SNP_select_file"
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
	if(sys.argv[2*i-1] == '--p_val'):
		p_val = float(sys.argv[2*i])
		label = 1
if(label == 0):
	print "Please provide the parameter --p_val"
	sys.exit()


remma_epi_select(plink_prefix, output_prefix, SNP_select_file,
num_test = num_test, p_value = p_val)
