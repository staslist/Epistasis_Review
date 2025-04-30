# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 13:57:26 2023

@author: staslist
"""
import logging
logging.basicConfig(level=logging.INFO)
import numpy as np
from gmat.gmatrix import agmat
from gmat.gmatrix import dgmat_as
from gmat.uvlmm.uvlmm_varcom import wemai_multi_gmat
from gmat.remma.remma_epiAA import remma_epiAA
from gmat.remma import annotation_snp_pos

# Step 1: Calculate the genomic relationship matrix
home_dir = './'
# home_dir = '/gpfs/group/home/slistopad/REMMA/data/epigen/'
# home_dir = '/gpfs/group/home/slistopad/REMMA/data/mouse/'
bed_file = home_dir + 'Pure_Multiplicative_TwoPairs_BaselineAlpha10_InteractionAlpha2_Chr1_CEU_SNP1000_IND1000_MAF005_02'  # the prefix for the plink binary file
# bed_file = home_dir + 'plink'
agmat(bed_file, home_dir + 'add_genom_rel_matrix_aa2_pure_mult_twopairs_interalpha2') 
dgmat_as(bed_file, home_dir + 'dom_genom_rel_matrix_aa2_pure_mult_twopairs_interalpha2')

# Step 2: Estimate the variances
add_genom_rel_matrix = home_dir + 'add_genom_rel_matrix_aa2_pure_mult_twopairs_interalpha2.agrm.mat_fmt'
dom_genom_rel_matrix = home_dir + 'dom_genom_rel_matrix_aa2_pure_mult_twopairs_interalpha2.dgrm_as.mat_fmt'
# pheno_file = home_dir + 'pheno_impute_median.ped'  # phenotypic file
pheno_file = home_dir + 'Pure_Multiplicative_TwoPairs_BaselineAlpha10_InteractionAlpha2_Chr1_CEU_SNP1000_IND1000_MAF005_02_Pheno.ped'
# pheno_file = home_dir + 'pheno'  # phenotypic file
ag = np.loadtxt(add_genom_rel_matrix)  # load the additive genomic relationship matrix
dg = np.loadtxt(dom_genom_rel_matrix)
#gmat_lst = [ag, dg, ag*ag, ag*dg, dg*dg]  # ag*ag is the additive by additive genomic relationship matrix
gmat_lst = [ag, dg, ag*ag]
wemai_multi_gmat(pheno_file, bed_file, gmat_lst, out_file= home_dir + 'var_aa2_pure_mult_twopairs_interalpha2.txt')


# Step 3: Test
var_com = np.loadtxt(home_dir + 'var_aa2_pure_mult_twopairs_interalpha2.txt')  # numpy array： [0] addtive variance; [1] additive by additive variance; [2] residual variance
remma_epiAA(pheno_file, bed_file, gmat_lst, var_com, p_cut=1.0e-5, out_file=home_dir + 'epiAA_aa2_pure_mult_twopairs_interalpha2')


# Step 4: Select top SNPs and add the SNP position
res_file = home_dir + 'epiAA_aa2_pure_mult_twopairs_interalpha2'  # result file
annotation_snp_pos(res_file, bed_file, p_cut=1.0e-5, dis=0)  # p values < 1.0e-5 and the distance between SNP pairs > 0
