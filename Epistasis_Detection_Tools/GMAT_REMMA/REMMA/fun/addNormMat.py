from pysnptools.snpreader import Bed
import numpy as np

def addNormMat(bed_file):
	
	snp_mat = Bed(bed_file, count_A1 = False).read()
	freq = np.sum(snp_mat.val, axis = 0) / (2*snp_mat.iid_count)
	freq.shape = (1, snp_mat.sid_count)
	snp_mat.val = np.subtract(snp_mat.val, 2*freq)
	return snp_mat.val.T

