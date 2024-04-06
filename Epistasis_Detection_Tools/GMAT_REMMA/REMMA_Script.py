# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 13:57:26 2023

@author: staslist
"""
import logging
logging.basicConfig(level=logging.INFO)
import numpy as np
from gmat.gmatrix import agmat
from gmat.uvlmm.uvlmm_varcom import wemai_multi_gmat
from gmat.remma.remma_epiAA import remma_epiAA
from gmat.remma.remma_epiAA import remma_epiAA_approx
from gmat.remma import annotation_snp_pos

# Step 1: Calculate the genomic relationship matrix
home_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch4_GWAS/Epistasis_Detection_Tools/GMAT_REMMA/gmat/data/native_american_DELETE/'
#home_dir = 'C:/Stas/LabWork/Bioinformatics/Projects/Ch4_GWAS/Epistasis_Detection_Tools/GMAT_REMMA/gmat/data/mouse/'
bed_file = home_dir + 'NA3_chr22'  # the prefix for the plink binary file
#bed_file = home_dir + 'plink'  # the prefix for the plink binary file
# agmat(bed_file, home_dir + 'add_genom_rel_matrix') 

# Step 2: Estimate the variances
add_genom_rel_matrix = home_dir + 'add_genom_rel_matrix.agrm.mat_fmt'
pheno_file = home_dir + 'pheno_impute_median.ped'  # phenotypic file
#pheno_file = home_dir + 'pheno'  # phenotypic file
ag = np.loadtxt(add_genom_rel_matrix)  # load the additive genomic relationship matrix
gmat_lst = [ag, ag*ag]  # ag*ag is the additive by additive genomic relationship matrix
# wemai_multi_gmat(pheno_file, bed_file, gmat_lst, out_file= home_dir + 'var_a_axa.txt')

# # Step 3: Test
var_com = np.loadtxt(home_dir + 'var_a_axa.txt')  # numpy array： [0] addtive variance; [1] additive by additive variance; [2] residual variance
#remma_epiAA(pheno_file, bed_file, gmat_lst, var_com, p_cut=1.0e-5, out_file=home_dir + 'epiAA_a_axa')
remma_epiAA_approx(pheno_file, bed_file, gmat_lst, var_com, p_cut=1.0e-5, out_file=home_dir + 'epiAA_approx_a_axa')

# # Step 4: Select top SNPs and add the SNP position
res_file = home_dir + 'epiAA_a_axa'  # result file
annotation_snp_pos(res_file, bed_file, p_cut=1.0e-5, dis=0)  # p values < 1.0e-5 and the distance between SNP pairs > 0
