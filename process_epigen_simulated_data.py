# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 04:31:19 2023

@author: staslist
"""

import numpy as np
import matplotlib.pyplot as plt

import inspect
import csv

import unittest

from statistics import mean

files = ['Impure_Multiplicative_OnePair_BaselineAlpha10_InteractionAlpha2_Chr1_CEU_SNP1000_IND1000_MAF005_02',
         'Impure_Multiplicative_OnePair_BaselineAlpha10_InteractionAlpha3_Chr1_CEU_SNP1000_IND1000_MAF005_02',
         'Impure_Multiplicative_OnePair_BaselineAlpha10_InteractionAlpha15_Chr1_CEU_SNP1000_IND1000_MAF005_02',
         'Impure_Multiplicative_OnePair_BaselineAlpha10_InteractionAlpha125_Chr1_CEU_SNP1000_IND1000_MAF005_02',
         'Impure_Multiplicative_TwoPairs_BaselineAlpha10_InteractionAlpha2_Chr1_CEU_SNP1000_IND1000_MAF005_02',
         'Impure_Multiplicative_TwoPairs_BaselineAlpha10_InteractionAlpha3_Chr1_CEU_SNP1000_IND1000_MAF005_02',
         'Impure_Multiplicative_TwoPairs_BaselineAlpha10_InteractionAlpha15_Chr1_CEU_SNP1000_IND1000_MAF005_02',
         'Impure_Multiplicative_TwoPairs_BaselineAlpha10_InteractionAlpha125_Chr1_CEU_SNP1000_IND1000_MAF005_02',
         'Pure_Multiplicative_OnePair_BaselineAlpha10_InteractionAlpha2_Chr1_CEU_SNP1000_IND1000_MAF005_02',
         'Pure_Multiplicative_OnePair_BaselineAlpha10_InteractionAlpha3_Chr1_CEU_SNP1000_IND1000_MAF005_02',
         'Pure_Multiplicative_OnePair_BaselineAlpha10_InteractionAlpha15_Chr1_CEU_SNP1000_IND1000_MAF005_02',
         'Pure_Multiplicative_OnePair_BaselineAlpha10_InteractionAlpha125_Chr1_CEU_SNP1000_IND1000_MAF005_02',
         'Pure_Multiplicative_TwoPairs_BaselineAlpha10_InteractionAlpha2_Chr1_CEU_SNP1000_IND1000_MAF005_02',
         'Pure_Multiplicative_TwoPairs_BaselineAlpha10_InteractionAlpha3_Chr1_CEU_SNP1000_IND1000_MAF005_02',
         'Pure_Multiplicative_TwoPairs_BaselineAlpha10_InteractionAlpha15_Chr1_CEU_SNP1000_IND1000_MAF005_02',
         'Pure_Multiplicative_TwoPairs_BaselineAlpha10_InteractionAlpha125_Chr1_CEU_SNP1000_IND1000_MAF005_02']

files_recessive = ['Impure_Recessive_OnePair_BaselineAlpha10_InteractionAlpha8_Chr1_CEU_SNP1000_IND1000_MAF025_04',
                   'Impure_Recessive_OnePair_BaselineAlpha10_InteractionAlpha16_Chr1_CEU_SNP1000_IND1000_MAF025_04',
                   'Impure_Recessive_TwoPairs_BaselineAlpha10_InteractionAlpha8_Chr1_CEU_SNP1000_IND1000_MAF025_04',
                   'Impure_Recessive_TwoPairs_BaselineAlpha10_InteractionAlpha16_Chr1_CEU_SNP1000_IND1000_MAF025_04',
                   'Pure_Recessive_OnePair_BaselineAlpha10_InteractionAlpha8_Chr1_CEU_SNP1000_IND1000_MAF025_04',
                   'Pure_Recessive_OnePair_BaselineAlpha10_InteractionAlpha16_Chr1_CEU_SNP1000_IND1000_MAF025_04',
                   'Pure_Recessive_TwoPairs_BaselineAlpha10_InteractionAlpha8_Chr1_CEU_SNP1000_IND1000_MAF025_04',
                   'Pure_Recessive_TwoPairs_BaselineAlpha10_InteractionAlpha16_Chr1_CEU_SNP1000_IND1000_MAF025_04']

files_dominant = ['Impure_Dominant_OnePair_BaselineAlpha10_InteractionAlpha8_Chr1_CEU_SNP1000_IND1000_MAF005_02',
                  'Impure_Dominant_OnePair_BaselineAlpha10_InteractionAlpha16_Chr1_CEU_SNP1000_IND1000_MAF005_02',
                  'Impure_Dominant_TwoPairs_BaselineAlpha10_InteractionAlpha8_Chr1_CEU_SNP1000_IND1000_MAF005_02',
                  'Impure_Dominant_TwoPairs_BaselineAlpha10_InteractionAlpha16_Chr1_CEU_SNP1000_IND1000_MAF005_02',
                  'Pure_Dominant_OnePair_BaselineAlpha10_InteractionAlpha8_Chr1_CEU_SNP1000_IND1000_MAF005_02',
                  'Pure_Dominant_OnePair_BaselineAlpha10_InteractionAlpha16_Chr1_CEU_SNP1000_IND1000_MAF005_02',
                  'Pure_Dominant_TwoPairs_BaselineAlpha10_InteractionAlpha8_Chr1_CEU_SNP1000_IND1000_MAF005_02',
                  'Pure_Dominant_TwoPairs_BaselineAlpha10_InteractionAlpha16_Chr1_CEU_SNP1000_IND1000_MAF005_02']

files_control = ['Pure_Control_TwoPairs_BaselineAlpha10_Chr1_CEU_SNP1000_IND1000_MAF005_02',
                 'Pure_Control_OnePair_BaselineAlpha10_Chr1_CEU_SNP1000_IND1000_MAF005_02',
                 'Pure_Control_EightPairs_BaselineAlpha10_Chr1_CEU_SNP1000_IND1000_MAF005_02',
                 'Impure_Control_TwoPairs_BaselineAlpha10_Chr1_CEU_SNP1000_IND1000_MAF005_02',
                 'Impure_Control_OnePair_BaselineAlpha10_Chr1_CEU_SNP1000_IND1000_MAF005_02',
                 'Impure_Control_EightPairs_BaselineAlpha10_Chr1_CEU_SNP1000_IND1000_MAF005_02']

files_bonus = ['Pure_XOR_EightPairs_BaselineAlpha10_InteractionAlpha16_Chr1_CEU_SNP1000_IND1000_MAF005_02',
               'Pure_XOR_EightPairs_BaselineAlpha10_InteractionAlpha8_Chr1_CEU_SNP1000_IND1000_MAF005_02',
               'Pure_Multiplicative_EightPairs_BaselineAlpha10_InteractionAlpha3_Chr1_CEU_SNP1000_IND1000_MAF005_02',
               'Pure_Multiplicative_EightPairs_BaselineAlpha10_InteractionAlpha2_Chr1_CEU_SNP1000_IND1000_MAF005_02',
               'Impure_XOR_EightPairs_BaselineAlpha10_InteractionAlpha16_Chr1_CEU_SNP1000_IND1000_MAF005_02',
               'Impure_XOR_EightPairs_BaselineAlpha10_InteractionAlpha8_Chr1_CEU_SNP1000_IND1000_MAF005_02',
               'Impure_Multiplicative_EightPairs_BaselineAlpha10_InteractionAlpha3_Chr1_CEU_SNP1000_IND1000_MAF005_02',
               'Impure_Multiplicative_EightPairs_BaselineAlpha10_InteractionAlpha2_Chr1_CEU_SNP1000_IND1000_MAF005_02']

def convert_raw_and_tfam_to_csv_and_pheno(raw_file:str, tfam_file:str):
    # Convert plink raw and fam files to csv and pheno files for use with Matrix Epistasis
    # Also generate the QMDR and MDR (binary phenotype) csv files
    # Also generate the epiSNP .txt and .dat files
    snp_names = []
    individuals_ordered = []
    csv_to_write = []
    with open(raw_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=' ')
        first_line = True
        for row in csv_reader:
            if(first_line):
                snp_names = row[6:]
                first_line = False
            else:
                individuals_ordered.append(row[0])
                csv_to_write.append(row[6:])
     
    pheno_dict = dict()
    pheno_dict_binary = dict()
    with open(tfam_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=' ')
        for row in csv_reader:
            pheno_dict[row[0]] = row[5]
            pheno = float(row[5])
            binary_pheno = 0
            if(pheno > 20.42):
                binary_pheno = 1
            pheno_dict_binary[row[0]] = str(binary_pheno)
            
    out_name = raw_file[0:-4]
    csv_name = out_name + '.csv'
    pheno_name = out_name + '.pheno'
    csv_name2 = out_name + '_QMDR.csv'
    csv_name3 = out_name + '_MDR.csv'
    txt_name = out_name + '_epiSNP.txt'
    dat_name = out_name + '_epiSNP.dat'
    epiSNP_id_map_file = out_name + '_epiSNP_ID_Map.csv'
    '''
    with open(csv_name, 'w') as writer:
        delimiter = ","
        writer.write(delimiter.join(snp_names) + '\n')
        for row in csv_to_write:
            my_string = delimiter.join(row)
            writer.write(my_string + '\n')
            
    with open(pheno_name, 'w') as writer:
        for indiv in individuals_ordered:
            writer.write(pheno_dict[indiv] + '\n')
            
    with open(csv_name2, 'w') as writer:
        delimiter = ","
        writer.write(delimiter.join(snp_names) + ',phenotype\n')
        i = 0
        for row in csv_to_write:
            my_string = delimiter.join(row)
            writer.write(my_string + ',')
            writer.write(pheno_dict[individuals_ordered[i]] + '\n')
            i += 1
    '''       
    with open(csv_name3, 'w') as writer:
        delimiter = ","
        writer.write(delimiter.join(snp_names) + ',phenotype\n')
        i = 0
        for row in csv_to_write:
            my_string = delimiter.join(row)
            writer.write(my_string + ',')
            writer.write(pheno_dict_binary[individuals_ordered[i]] + '\n')
            i += 1
    '''
    epiSNP_id_map = dict()
    with open(dat_name, 'w') as writer:
        writer.write('fam_ID\tSex\tind_ID\t')
        delimiter = ","
        writer.write(delimiter.join(snp_names) + '\n')
        i = 0
        for row in csv_to_write:
            epiSNP_id_map[i] = individuals_ordered[i]
            #writer.write(individuals_ordered[i] + '\t9\t' + individuals_ordered[i] + '\t')
            writer.write(str(i) + '\t-1\t' + str(i) + '\t')
            my_string = delimiter.join(row)
            writer.write(my_string + '\n')
            i += 1
    
    with open(epiSNP_id_map_file, 'w') as writer:
        for k,v in epiSNP_id_map.items():
            writer.write(str(k) + ',' + v + '\n')
    
    with open(txt_name, 'w') as writer:
        writer.write('ind_ID\tfather_ID\tmother_ID\tSex\tPhenotype\n')
        i = 0
        for indiv in individuals_ordered:
            father_id = str(len(individuals_ordered) + i)
            mother_id = str(len(individuals_ordered)*2 + i)
            writer.write(str(i) + '\t' + father_id + '\t' + mother_id + '\t-1\t' + pheno_dict[indiv] + '\n')
            i += 1
    '''
def convert_fam_to_pheno(curr_dir:str, fam_file:str):
    to_write = []
    with open(curr_dir + '/' + fam_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            to_write.append(row)
            
    fname_out = fam_file[:-4] + '_Pheno.ped'
    with open(curr_dir + '/' + fname_out, 'w') as f:
        f.write('famid' + '\t' + 'id' + '\t' + 'mean' + '\t' + 'trait\n')
        for row in to_write:
            line = row[0] + '\t' + row[1] + '\t1\t' + row[5] + '\n' 
            f.write(line)
            

def check_two_files_identical(file1:str, file2:str):
    lines1 = []
    with open(file1, 'r') as f:
        lines1 = f.readlines()
        
    lines2 = []
    with open(file2, 'r') as f:
        lines2 = f.readlines()
        
    i = 0
    for line in lines1:
        if(line != lines2[i]):
            return False
        i += 1
    return True
    

def check_epigen_json_for_detectability(fname_json:str, inter_dis_snp_indeces:list, interaction:str):
    # This function checks whether the recessive interaction is possible to detect.
    # That is there must be at least some individual for whom two interacting disease snps both have genotype 2.
    
    num_individuals = 0
    num_snps = 0
    phenotypes = []
    individual_genotypes = []
    individual_snps = []
    
    disease_snp_names = []
    interacting_disease_snp_names = []
    interacting_disease_snp_indeces = []
    
    with open(fname, 'r') as f:
        lines = f.readlines()
        #print(lines)
        
        #print(len(lines[0]))
        num_snps_index = lines[0].find('"num_snps"')
        num_inds_index = lines[0].find('"num_inds"')
        model_type_index = lines[0].find('"model_type"')
        
        # NUMBER OF INDIVIDUALS PROCESSING
        
        num_individuals = int(lines[0][num_inds_index+12:model_type_index-2])
        
        genotype_index = lines[0].find('"genotype"')
        phen_index = lines[0].find('"phenotype"')
        snps_index = lines[0].find('"snps"')
        disease_snps_index = lines[0].find('"disease_snps"')
        mafs_index = lines[0].find('"mafs"')
        
        # GENOTYPE PROCESSING
        genotype = lines[0][genotype_index+13:phen_index-3]
        temp = ''
        append_flag = False
        #print(genotype)
        for char in genotype:
            if(char == '['):
                append_flag = True
                continue
            elif(char == ']'):
                append_flag = False
                individual_genotypes.append(temp.split(', '))
                temp = ''
                continue
            if(append_flag):
                temp = temp + char
        
        # SNP DATA PROCESSING
        snps = lines[0][snps_index + 9:disease_snps_index-3]
        
        temp = ''
        append_flag = False
        #print(phenotypes_array)
        for char in snps:
            if(char == '['):
                append_flag = True
                continue
            elif(char == ']'):
                append_flag = False
                num_snps += 1
                individual_snps.append(temp.replace('"', '').split(', '))
                temp = ''
                continue
            if(append_flag):
                temp = temp + char
        
        # print(individual_snps)
        # REPORT DISEASE SNPS
        disease_snps = lines[0][disease_snps_index+17:mafs_index-3]
        disease_snps = disease_snps.split(', ')
        #print("DISEASE SNPs in file: ", fname)
        #print(disease_snps)
        for ds_snp in disease_snps:
            #print(individual_snps[int(ds_snp)])
            disease_snp_names.append(individual_snps[int(ds_snp)])
            
        for tuple_pair in inter_dis_snp_indeces:
            interacting_disease_snp_names.append((disease_snp_names[tuple_pair[0]][0], disease_snp_names[tuple_pair[1]][0]))
            interacting_disease_snp_indeces.append((int(disease_snps[tuple_pair[0]]), int(disease_snps[tuple_pair[1]])))
        
        #print(interacting_disease_snp_indeces)
        
        num_pairs = len(interacting_disease_snp_indeces)
        
        detectable_pairs = 0
        # detectable = False
        index = 0
        for inter_dis_snp_pair in interacting_disease_snp_indeces:
            i = 0
            #print('Pair ', str(index))
            #print(individual_genotypes[inter_dis_snp_pair[0]])
            #print(individual_genotypes[inter_dis_snp_pair[1]])
            while i < num_individuals:
                if(interaction == 'Recessive'):
                    if(individual_genotypes[inter_dis_snp_pair[0]][i] == '2' and 
                       individual_genotypes[inter_dis_snp_pair[1]][i] == '2'):
                        detectable_pairs += 1
                        break
                elif(interaction == 'Dominant'):
                    if(individual_genotypes[inter_dis_snp_pair[0]][i] in ['1','2'] and 
                       individual_genotypes[inter_dis_snp_pair[1]][i] in ['1','2']):
                        detectable_pairs += 1
                        break
                elif(interaction == 'Multiplicative'):
                    if(individual_genotypes[inter_dis_snp_pair[0]][i] in ['1', '2'] or 
                       individual_genotypes[inter_dis_snp_pair[1]][i] in ['1','2']):
                        detectable_pairs += 1
                        break
                elif(interaction == 'XOR'):
                    if((individual_genotypes[inter_dis_snp_pair[0]][i] in ['1', '2'] and 
                       individual_genotypes[inter_dis_snp_pair[1]][i] == '0') or 
                       (individual_genotypes[inter_dis_snp_pair[1]][i] in ['1', '2'] and 
                        individual_genotypes[inter_dis_snp_pair[0]][i] == '0')):
                        detectable_pairs += 1
                        break
                    
                i += 1
            index += 1
            
        print("Epigen JSON filename: ", fname_json)
        print("Was " + interaction + " interaction detectable:", bool(detectable_pairs))
        print("Number of disease snp pairs: ", num_pairs)
        print("Number of detectable interacting pairs: ",  detectable_pairs)
   
def convert_epigent_qt_to_cc(fname:str):
    ''' Convert input simulated epigen data with quantitative phenotype to case control phenotype.
    Any phenotype with value > threshold. Any phenotype with value < threshold is control.
    The threshold is custom.
    Assume all simulated data exists as one line (first line in the file).'''
    
    threshold = 1
    
    if ('Multiplicative' in fname):
        if('InteractionAlpha125' in fname):
            threshold = 15
        elif('InteractionAlpha15' in fname):
            threshold = 22
        elif('InteractionAlpha2' in fname):
            threshold = 40
        elif('InteractionAlpha3' in fname):
            threshold = 90
        else:
            raise ValueError('Invalid interaction alpha.')
    elif('Control' in fname):
        threshold = 20
    else:
        if('InteractionAlpha8' in fname):
            threshold = 30
        elif('InteractionAlpha16' in fname):
            threshold = 60
        else:
            raise ValueError('Invalid interaction alpha.')
    
    output = ''
    with open(fname, 'r') as f:
        lines = f.readlines()
        # print(lines)
        
        # print(len(lines[0]))
        num_snps_index = lines[0].find('"num_snps"')
        num_inds_index = lines[0].find('"num_inds"')
        model_type_index = lines[0].find('"model_type"')
        genotype_index = lines[0].find('"genotype"')
        phen_index = lines[0].find('"phenotype"')
        snps_index = lines[0].find('"snps"')
        disease_snps_index = lines[0].find('"disease_snps"')
        mafs_index = lines[0].find('"mafs"')
        

        # PHENOTYPE PROCESSING
        pre_phe_content = lines[0][0:phen_index+13]
        post_phen_content = lines[0][snps_index-2:]
        phenotype = lines[0][phen_index+14:snps_index-3]
        phenotypes = phenotype.split(', ')
        
        # print(phenotypes)
        
        phenotypes_float = []
        for ph in phenotypes:
            phenotypes_float.append(float(ph))
        
        # phenotypes_array = np.array(phenotypes_float)
        # pheno_max = np.max(phenotypes_array)
        # pheno_min = np.min(phenotypes_array)
        # threshold = (pheno_max + pheno_min) * factor
        
        phenotypes_converted = []
        # Control: 0 Case: 1
        case_counter = 0
        for pheno in phenotypes_float:
            if(pheno > threshold):
                phenotypes_converted.append(2)
                case_counter += 1
            else:
                phenotypes_converted.append(1)
                
        # print(case_counter)
                
        output = pre_phe_content + str(phenotypes_converted) + post_phen_content
        # print(output)
        
    fname_list = fname.split('.')
    fname_out = fname_list[0] + '_converted.json'
    
    with open(fname_out, 'w') as f2:
        f2.write(output) 
   
# def convert_epigens_to_ped_and_map(fnames:list): 
#     ''' One time use function to generate extra large .ped and .map files. '''
   
#     num_individuals, num_snps, phenotypes, individual_genotypes, individual_snps = extract_key_info_from_multiple_EpiGEN_jsons(fnames)
    
#     #print(phenotypes)
#     #print(individual_genotypes)
#     #print(individual_snps)
    
#     fname_out_ped = 'Chrs1_to_12_CEU_SNP1000000_IND1000_MAF005_02.ped'
#     fname_out_map = 'Chrs1_to_12_CEU_SNP1000000_IND1000_MAF005_02.map'
    
#     with open(fname_out_ped, 'w') as f2:
#         f2.write('#Family_ID	Individual_ID	Paternal_ID	Maternal_ID	Sex	Phenotype	SNP1	SNP2	...	SNPN\n')
#         i = 1
#         while i <= num_individuals:
#             f2.write('FAM' + str(i) + '\t' + str(i) + '\t' + '0\t0\t0\t' + str(phenotypes[i-1]) + '\t')
#             j = 0
#             while j < num_snps:
#                 geno = individual_genotypes[j][i-1]
#                 major_allele = individual_snps[j][3]
#                 minor_allele = individual_snps[j][4]
#                 if(geno == '0'):
#                     f2.write(str(major_allele) + '\t' + str(major_allele))
#                 elif(geno == '1'):
#                     f2.write(str(major_allele) + '\t' + str(minor_allele))
#                 elif(geno == '2'):
#                     f2.write(str(minor_allele) + '\t' + str(minor_allele))
#                 else:
#                     raise ValueError("Invalid genotype value in an epigen simulated data file.")
                    
#                 if(j < (num_snps-1)):
#                     f2.write('\t')
#                 j += 1
#             i += 1
#             f2.write('\n')
            
#     with open(fname_out_map, 'w') as f3:
#         f3.write('#chr	snp	pos\n')
#         i = 0
#         while i < num_snps:
#             f3.write(individual_snps[i][1][3:] + '\t' + individual_snps[i][0] + '\t' + individual_snps[i][2] + '\n')
#             i += 1
    
def convert_epigen_to_DMM(fname:str, name:str):
    ''' This converts epigen simulation data (in json format) into an X and Y .npy matrices
    that can be used by DMM tool for marginal epistasis analysis. 
    X is a n x p matrix, wherein n is samples, and p is snps.
    Y is a n x 1 matrix, wherein in is samples, one phenotype is present.'''
    num_individuals, num_snps, phenotypes, genotypes, individual_snps = extract_key_info_from_EpiGEN_json(fname)
    #print(genotypes)
    #print(phenotypes)
    #print(individual_snps)
    
    # The Y matrix is already prepared (phenotypes variable)
    # The individual genotypes is currently an array of arrays, wherein each subarray is 
    # genotypes of all individuals for one SNP. This needs to be re-arranged into an X matrix.
    
    X = []
    j = 0
    limit = len(genotypes[0])
    while j < limit:
        row = []
        for geno in genotypes:
            row.append(geno[j])
            
        X.append(row)
        j += 1
        
    print(X)
    X = np.array(X)
    
    path = '...'
    pathX = path + name + '_X'
    pathY = path + name + '_Y'
    np.save(pathX, X)
    np.save(pathY, phenotypes)
    

def convert_epigen_to_ped_and_map(fname:str):
    ''' This converts epigen simulation data (in json format) to .ped and .map files necessary for plink 1.9 input.
    If write is false this function instead returns the disease and interacting disease snps.'''
    # inter_dis_snp_indeces specifies the indeces of interacting disease snps within the disease snp list.
    # currently assumes that all interacting snps are pairs.
    
    num_individuals, num_snps, phenotypes, individual_genotypes, individual_snps = extract_key_info_from_EpiGEN_json(fname)
    
    #print(phenotypes)
    #print(individual_genotypes)
    #print(individual_snps)
    
    fname_list = fname.split('.')
    fname_out_ped = fname_list[0] + '.ped'
    fname_out_map = fname_list[0] + '.map'
    
    with open(fname_out_ped, 'w') as f2:
        f2.write('#Family_ID	Individual_ID	Paternal_ID	Maternal_ID	Sex	Phenotype	SNP1	SNP2	...	SNPN\n')
        i = 1
        while i <= num_individuals:
            f2.write('FAM' + str(i) + '\t' + str(i) + '\t' + '0\t0\t0\t' + str(phenotypes[i-1]) + '\t')
            j = 0
            while j < num_snps:
                geno = individual_genotypes[j][i-1]
                major_allele = individual_snps[j][3]
                minor_allele = individual_snps[j][4]
                if(geno == '0'):
                    f2.write(str(major_allele) + '\t' + str(major_allele))
                elif(geno == '1'):
                    f2.write(str(major_allele) + '\t' + str(minor_allele))
                elif(geno == '2'):
                    f2.write(str(minor_allele) + '\t' + str(minor_allele))
                else:
                    raise ValueError("Invalid genotype value in an epigen simulated data file.")
                    
                if(j < (num_snps-1)):
                    f2.write('\t')
                j += 1
            i += 1
            f2.write('\n')
            
    with open(fname_out_map, 'w') as f3:
        f3.write('#chr	snp	pos\n')
        i = 0
        while i < num_snps:
            f3.write(individual_snps[i][1][3:] + '\t' + individual_snps[i][0] + '\t' + individual_snps[i][2] + '\n')
            i += 1
             
def convert_ped_and_map_to_SNPassoc_csv(fname_ped, fname_map):
    snp_names = []
    with open(fname_map, 'r') as f:
        lines = f.readlines()
        
        for line in lines:
            split_line = line.split('\t')
            snp_names.append(split_line[1])
    
    genotype_info = []
    with open(fname_ped, 'r') as f:
        lines = f.readlines()
        
        first_line = True
        for line in lines:
            if(not first_line):
                split_line = line.split('\t')
                temp_start = [split_line[0], split_line[5]]
                genotypes = split_line[6:]
                #print(len(genotypes))
                genotypes_paired = []
                i = 0
                while i < len(genotypes):
                    genotypes_paired.append(genotypes[i][0] + '' + genotypes[i+1][0])
                    i += 2
                genotype_info.append(temp_start + genotypes_paired)
            
            first_line = False
    #print(genotype_info)
    
    header_line = ['ID', 'Phenotype'] + snp_names[1:]
    #print(header_line)
    #print(genotype_info)
    
    with open(fname_ped[0:-4] + '_SNPassoc.csv', 'w') as f:
        i = 0
        for ele in header_line:
            if(i == (len(header_line)-1)):
                f.write(ele)
            else:
                f.write(ele + '\t')
            i += 1
        f.write('\n')
        
        z = 0
        for row in genotype_info:
            i = 0
            for ele in row:
                if(i == (len(row) - 1)):
                    f.write(ele)
                else:
                    f.write(ele + '\t')
                i += 1
            if(z != (len(genotype_info) - 1)):
               f.write('\n')
            z += 1

def extract_key_info_from_multiple_EpiGEN_jsons(fname_jsons:list):
    num_inds_total = 0
    num_snps_total = 0
    phenotypes_total = np.array([])
    ind_genotypes_total = []
    ind_snps_total = []
    for fname_json in fname_jsons:
        num_inds, num_snps, phenotypes, ind_genotypes, ind_snps = extract_key_info_from_EpiGEN_json(fname_json)
        num_inds_total += num_inds
        num_snps_total += num_snps
        phenotypes_total = np.concatenate((phenotypes_total, phenotypes))
        for ind_genotype in ind_genotypes:
            ind_genotypes_total.append(ind_genotype)
        for ind_snp in ind_snps:
            ind_snps_total.append(ind_snp)
            
    return num_inds_total, num_snps_total, phenotypes_total, ind_genotypes_total, ind_snps_total

def extract_key_info_from_EpiGEN_json(fname_json:str):
    ''' Return  '''
    
    num_individuals = 0
    num_snps = 0
    phenotypes = []
    individual_genotypes = []
    individual_snps = []
    
    with open(fname_json, 'r') as f:
        lines = f.readlines()
        #print(lines)
        
        #print(len(lines[0]))
        num_snps_index = lines[0].find('"num_snps"')
        num_inds_index = lines[0].find('"num_inds"')
        model_type_index = lines[0].find('"model_type"')
        
        # NUMBER OF INDIVIDUALS PROCESSING
        
        num_individuals = int(lines[0][num_inds_index+12:model_type_index-2])
        #print(num_individuals)
        
        genotype_index = lines[0].find('"genotype"')
        phen_index = lines[0].find('"phenotype"')
        snps_index = lines[0].find('"snps"')
        disease_snps_index = lines[0].find('"disease_snps"')

        # PHENOTYPE PROCESSING
        phenotype = lines[0][phen_index+14:snps_index-3]
        phenotypes = phenotype.split(', ')
        phenotypes_float = []
        for ph in phenotypes:
            phenotypes_float.append(float(ph))
        
        phenotypes_array = np.array(phenotypes_float)
        # print("Max:")
        # print(np.max(phenotypes_array))
        # print("Min:")
        # print(np.min(phenotypes_array))
        # print("Mean:")
        # print(np.mean(phenotypes_array))
        # print("STD:")
        # print(np.std(phenotypes_array))
        
        #plt.hist(phenotypes_array, range = (0,100))
        # plt.hist(phenotypes_array)
        
        # GENOTYPE PROCESSING
        genotype = lines[0][genotype_index+13:phen_index-3]
        temp = ''
        append_flag = False
        #print(genotype)
        for char in genotype:
            if(char == '['):
                append_flag = True
                continue
            elif(char == ']'):
                append_flag = False
                individual_genotypes.append(temp.split(', '))
                temp = ''
                continue
            if(append_flag):
                temp = temp + char
        
        #print(individual_genotypes)
        
        # SNP DATA PROCESSING
        snps = lines[0][snps_index + 9:disease_snps_index-3]
        
        temp = ''
        append_flag = False
        #print(phenotypes_array)
        for char in snps:
            if(char == '['):
                append_flag = True
                continue
            elif(char == ']'):
                append_flag = False
                num_snps += 1
                individual_snps.append(temp.replace('"', '').split(', '))
                temp = ''
                continue
            if(append_flag):
                temp = temp + char
        
        #print(num_snps)
        #print(individual_snps)    

    return num_individuals, num_snps, phenotypes_array, individual_genotypes, individual_snps

def convert_epigen_json_to_EpiSNP_files(fname_json:str):
    ''' Need to create two files. chrxxx.dat and insert_name.txt 
    The former file contains following columns: fam_ID, Sex, ind_ID, snp1, snp2, etc...
    The latter file contains following columns: ind_ID, father_ID, mother_ID, Sex, Phenotype.
    Assume # of individuals = N.
    Assume all individuals are unrelated. Have individual and family IDs go from 1 to 
    N. Have father ids go from (N + 1) to 2 * (N). Have mother ids go from (2*N + 1) to (3 * N).
    Have sex for all individuals = -1 (AKA unkonwn).
    Assume that there is only one chromosome file (presumably chromosome 1).'''
    
    num_samples, num_snps, phenotypes, genotypes, snps = extract_key_info_from_EpiGEN_json(fname_json)
    
    out_fname1 = fname_json[0:-5] + '_chr1.dat'
    out_fname2 = fname_json[0:-5] + '.txt'
    
    short_name = ''
    if('Impure' in fname_json):
        short_name += 'I'
    else:
        short_name += 'P'
        
    if('Dominant' in fname_json):
        short_name += '_D'
    elif('Multiplicative' in fname_json):
        short_name += '_M'
    elif('Recessive' in fname_json):
        short_name += '_R'
    elif('XOR' in fname_json):
        short_name += '_X'
        
    if('OnePair' in fname_json):
        short_name += '_OP'
    elif('TwoPairs' in fname_json):
        short_name += '_TP'
    elif('EightPairs' in fname_json):
        short_name += '_EP'
        
    if('InteractionAlpha125' in fname_json):
        short_name += '_IA125'
    elif('InteractionAlpha15' in fname_json):
        short_name += '_IA15'
    elif('InteractionAlpha2' in fname_json):
        short_name += '_IA2'
    elif('InteractionAlpha3' in fname_json):
        short_name += '_IA3'
    elif('InteractionAlpha8' in fname_json):
        short_name += '_IA8'
    elif('InteractionAlpha16' in fname_json):
        short_name += '_IA16'
        
    if('MAF005_02' in fname_json):
        short_name += '_MAF005_02'
    elif('MAF025_04' in fname_json):
        short_name += '_MAF025_04'
        
    short_name1 = short_name + '_chr1.dat'
    short_name2 = short_name + '.txt'
    
    
    with open(short_name1, 'w') as f:
        i = 1
        f.write('fam_ID\tSex\tind_ID\t')
        for snp in snps:
            f.write(snp[0])
            if(i < len(snps)):
                f.write(',')
                
            i += 1
        f.write('\n')
        
        j = 0
        limit = len(genotypes[0])
        while j < limit:
            # Write fam_ID, Sex, and ind_ID values
            f.write(str(j) + '\t-1\t' + str(j) + '\t')
            k = 1
            for geno in genotypes:
                f.write(geno[j])
                if(k < len(genotypes)):
                    f.write(',')
                k += 1
                
            f.write('\n')
            j += 1
            
    with open(short_name2, 'w') as f:
        f.write('ind_ID\tfather_ID\tmother_ID\tSex\tPhenotype\n')
        j = 0
        for phenotype in phenotypes:
            f.write(str(j) + '\t' + str(num_samples + j) + '\t' + str(2*num_samples + j) + '\t-1\t')
            f.write(str(phenotype))
            f.write('\n')
            j += 1

def convert_epigen_json_to_genotype_csv(fname_json:str, single_file:bool = False):
    ''' This converts epigen simulation data (in json format) to csv format for use with MatrixEpistasis 
    or QMDR.'''
    # inter_dis_snp_indeces specifies the indeces of interacting disease snps within the disease snp list.
    # currently assumes that all interacting snps are pairs.
    
    num_individuals, num_snps, phenotypes, individual_genotypes, individual_snps = extract_key_info_from_EpiGEN_json(fname_json)
    
    
    if(not single_file):
        out_fname = fname_json[0:-5] + '.csv'
    else:
        out_fname = fname_json[0:-5] + '_QMDR.csv'
    
    with open(out_fname, 'w') as f:
        i = 1
        for snp in individual_snps:
            f.write(snp[0])
            if(i < len(individual_snps)):
                f.write(',')
            elif(i == len(individual_snps) and single_file):
                f.write(',phenotype')
            
            i += 1
        f.write('\n')
        
        j = 0
        limit = len(individual_genotypes[0])
        while j < limit:
            k = 1
            for geno in individual_genotypes:
                f.write(geno[j])
                if(k < len(individual_genotypes)):
                    f.write(',')
                elif(k == len(individual_genotypes) and single_file):
                    f.write(',')
                    f.write(str(phenotypes[j]))
                k += 1
                
            f.write('\n')
            j += 1
    
    if(not single_file):
        with open(fname_json[0:-5] + '.pheno', 'w') as f:
            for phenotype in phenotypes:
                f.write(str(phenotype))
                f.write('\n')

def extract_interaction_pos(file_path):
	""" Get list of positions of interacting SNP groups as tuples from input model .xml file. """

	# Read file 
	with open(file_path, 'r') as file:
		lines = file.readlines()
	interaction_pos_list = []
	within_interaction_model = False
	current_pos_list = []
	# Iterate throughlines
	for line in lines:
		line = line.strip()
		if '<InteractionModel' in line:
			# New interacting SNP group starts
			within_interaction_model = True
		elif '</InteractionModel>' in line:
			# End of interacting SNP group
			if current_pos_list:
				# Add tuple of group list to overall list 
				interaction_pos_list.append(tuple(current_pos_list))
			current_pos_list = []
			within_interaction_model = False
		elif within_interaction_model and line.startswith('<pos>') and line.endswith('</pos>'):
			# Add position of SNP in group to group list
			pos_value = int(line[5:-6])
			current_pos_list.append(pos_value)
	return interaction_pos_list
            
def report_epistasis_results_helper2(sorted_epi_results: list, inter_dis_snps:set, num_tps:int, metric:str,
									 verbose:bool = True):
    """Returns true positive count, false positive count, false positive count, and penalized 
    and unpenalized average true positive position given results of tool sorted by metric, 
    list of interacting snps, and number of true positives.
	Code loops through sorted results and checks against list of interacting SNP pairs. 
	Prints results if verbose is True. 
	Code is modified from my mentor's code for reporting results from REMMA.
	"""
    true_pos_count = 0
    false_pos_count = 0
    index = 1
    tp_indeces = []
    tp_indeces_penalized = []
	# added set of true positives to make sure duplicates are not recorded
    tp_pairs = set()
    for tup in sorted_epi_results:
        match_flag = False
	
        # iterate through interacting pairs *not detected* (modification) 
        for pair in (inter_dis_snps-tp_pairs):
            rev_pair = (pair[1], pair[0])
            if(tup[0] == pair or tup[0] == rev_pair):
                true_pos_count += 1
                tp_indeces.append(index)
                tp_indeces_penalized.append(index)
                match_flag = True
				# add pair to set 
                tp_pairs.add(pair)
		
        if(not match_flag):
            false_pos_count += 1
        index += 1
	
    if(len(tp_indeces) == 0):
        tp_indeces = [0]
        tp_indeces_penalized = [0]
	# This elif penalizes missing true positives for datasets with more than one interacting snp pairs, 
	# by listing the undetected true positives as having the last position in the sorted list (mentor's documentation) 
    elif(len(tp_indeces) > 0 and len(tp_indeces) < num_tps):
        while(len(tp_indeces_penalized) < num_tps):
            tp_indeces_penalized.append(index)
		
    false_neg_count = num_tps - true_pos_count
    f1 = (2 * true_pos_count) / (2 * true_pos_count + false_pos_count + false_neg_count)
    
    if(verbose):
		#print(tp_pairs)
        print("Number of true positives: ", true_pos_count)
        print("Number of false positives: ", false_pos_count)
        print("F1 Score: ", f1)
        print("Average true positive position (unpenalized) ranked by " + metric + ": ", mean(tp_indeces))
        print("Average true positive position (penalized) ranked by " + metric + ": ", mean(tp_indeces_penalized))
	
    return true_pos_count, false_pos_count, mean(tp_indeces), mean(tp_indeces_penalized)
 
def report_remma_epistasis_results(fname_json:str, fname_xml:str, sim_data:str, return_dict = False, print_results=True):
    """ Given dataset information (json, model, and folder identifier), report results for all five configurations. 
    Return dictionary if return_dict is True, else return statistics. Print results if print_results is True. 
	"""
    # Get positions of interacting SNP pairs
    inter_dis_snp_indeces = extract_interaction_pos(fname_xml)
	
    #print(inter_dis_snp_indeces)
	# Get SNP numbers from positions (mentor's function)
    inter_dis_snps = get_disease_snps_from_epigen_json(fname_json, inter_dis_snp_indeces)
    #print(inter_dis_snps)
	# Number of true positives is number of pairs
    num_tps = len(inter_dis_snps)
	#print(inter_dis_snps)

    config_list = ['aa1', 'aa2', 'aa3', 'ad', 'dd']
    epi_results = {}
	# Iterate through all configurations 
    base_path = '...'
    for config in config_list: 
		# Get path to output file 
        if config == 'aa1':
            path = base_path + config+'/'+config+'_'+sim_data+'/epiAA_'+config+'_'+sim_data+'.anno' 
        elif config == 'aa2' or config == 'aa3':
            path = base_path + config+'/'+config+'_'+sim_data+'/epiAA_'+config+'_'+sim_data+'.anno'
        else:
            path = base_path + config+'/'+config+'_'+sim_data+'/epi'+config.upper()+'_'+config+'_'+sim_data+'.anno' 
		# Read file contents and get results
        with open(path, 'r') as f:
            lines = f.readlines()
			# Loop through lines
            for i in range(1, len(lines)):
                line = lines[i].split(' ')
                pval = float(line[-1][:-1])
                if (pval >= 1e-7):
                    continue
				# Extract corresponding pair
                key = (line[2], line[9])
				#print(key)
				# Add pair and p-value to results
                if key in epi_results:
					# Handle overlapping results between configurations by recording lowest p-value
                    epi_results[key] = min(epi_results[key], float(line[-1][:-1]))
                else:
                    epi_results[key] = float(line[-1][:-1])
	# Sort results by p-value
    sorted_epi_results = sorted(epi_results.items(), key=lambda x:x[1])

    #print(sorted_epi_results)
    if (return_dict):
		# Return sorted results
        return sorted_epi_results
    else: 
		# Calculate and return statistics 
        a, b, c, d = report_epistasis_results_helper2(sorted_epi_results, set(inter_dis_snps), num_tps, 'p-value', print_results)
        return (a, len(inter_dis_snps), b, c, d)
           
def report_epistasis_results_helper(sorted_epi_results: list, inter_dis_snps:list, num_tps:int, metric:str,
                                    verbose:bool = True):   
    true_pos_count = 0
    false_pos_count = 0
    index = 1
    tp_indeces = []
    tp_indeces_penalized = []
    #print(inter_dis_snps)
    for tup in sorted_epi_results:
        match_flag = False
        
        #print(tup)
        
        for pair in inter_dis_snps:
            rev_pair = (pair[1], pair[0])
            #print(tup[0])
            #print(pair)
            #print(rev_pair)
            if(tup[0] == pair or tup[0] == rev_pair):
                true_pos_count += 1
                tp_indeces.append(index)
                tp_indeces_penalized.append(index)
                #print(tup[0])
                #print(tup[1])
                #print(index)
                match_flag = True
        
        if(not match_flag):
            false_pos_count += 1
            
        index += 1
    
    #print(tp_indeces)
    if(len(tp_indeces) == 0):
        tp_indeces = [0]
        tp_indeces_penalized = [0]
    # This elif penalizes missing true positives for datasets with more than one interacting snp pairs, 
    # by listing the undetected true positives as having the last position in the sorted list. 
    elif(len(tp_indeces) > 0 and len(tp_indeces) < num_tps):
        while(len(tp_indeces_penalized) < num_tps):
            tp_indeces_penalized.append(index)
     
    false_neg_count = num_tps - true_pos_count
    f1 = (2 * true_pos_count) / (2 * true_pos_count + false_pos_count + false_neg_count)    
     
    if(verbose):
        print("Number of true positives: ", true_pos_count)
        print("Number of false positives: ", false_pos_count)
        print("F1 Score: ", f1)
        print("Average true positive position (unpenalized) ranked by " + metric + ": ", mean(tp_indeces))
        print("Average true positive position (penalized) ranked by " + metric + ": ", mean(tp_indeces_penalized))
    
    return true_pos_count, false_pos_count, mean(tp_indeces), mean(tp_indeces_penalized)
    
def report_MIDESP_epistasis_results(fname_json:str, inter_dis_snp_indeces:list, fname_epistasis:str,
                                    return_dict = False):
    inter_dis_snps = get_disease_snps_from_epigen_json(fname_json, inter_dis_snp_indeces)
    num_tps = len(inter_dis_snps)
    #print(inter_dis_snps)
    epi_results_dict = {}
    sorted_epi_results = {}
    with open(fname_epistasis, 'r') as f:
        lines = f.readlines()
        i = 1
        while i < len(lines):
            line = lines[i].split(' ')

            if(not (line[2], line[0]) in epi_results_dict):
                epi_results_dict[(line[0], line[2])] = float(line[3])
            i += 1
        sorted_epi_results = sorted(epi_results_dict.items(), key=lambda x:x[1], reverse=True)
        #print(sorted_epi_results)
    
    if(return_dict):
        return sorted_epi_results
    else:
        report_epistasis_results_helper(sorted_epi_results, inter_dis_snps, num_tps, 'MI')
  
def report_boost_epistasis_results(fname_json:str, inter_dis_snp_indeces:list, fname_epistasis:str, 
                                   return_dict:bool = False):  
    inter_dis_snps = get_disease_snps_from_epigen_json(fname_json, inter_dis_snp_indeces)
    num_tps = len(inter_dis_snps)
    # print(inter_dis_snps)
    
    epi_results_dict = {}
    sorted_epi_results = {}
    with open(fname_epistasis, 'r') as f:
        lines = f.readlines()
        i = 1
        while i < len(lines):
            #print(lines[i])
            
            line = lines[i].split(' ')
            line_cleaned = []
            for ele in line:
                if(ele != '' and ele != '\n'):
                    line_cleaned.append(ele)
            
            #print(line_cleaned)

            if(not (line_cleaned[3], line_cleaned[1]) in epi_results_dict):
                epi_results_dict[(line_cleaned[1], line_cleaned[3])] = float(line_cleaned[6])
            i += 1
        sorted_epi_results = sorted(epi_results_dict.items(), key=lambda x:x[1])
        #print(sorted_epi_results)
        
    if(return_dict):
        return sorted_epi_results
    else:
        report_epistasis_results_helper(sorted_epi_results, inter_dis_snps, num_tps, 'p-value')

def report_qmdr_epistasis_results(fname_json:str, inter_dis_snp_indeces:list, fname_epistasis:str,
                                  return_dict:bool = False):  
    inter_dis_snps = get_disease_snps_from_epigen_json(fname_json, inter_dis_snp_indeces)
    num_tps = len(inter_dis_snps)
    # print(inter_dis_snps)
    
    epi_results_dict = {}
    sorted_epi_results = {}
    key_lines = []
    with open(fname_epistasis, 'r') as f:
        lines = f.readlines()
        i = 1
        record_lines = False
        while i < len(lines):
            if(lines[i][0:8] == '@Results'):
                if(lines[i+1][0:14] == 'BeginTopModels'):
                    i += 2
                    record_lines = True
            elif(lines[i][0:12] == 'EndTopModels'):
                record_lines = False
                
            if(record_lines):
                key_lines.append(lines[i])

            i += 1
            
    for line in key_lines:
        snp_one = line.split(',')[0]
        part_two = line.split(',')[1]
        snp_two = part_two.split('\t')[0]
        bal_acc = part_two.split('\t')[1]

        if(not (snp_one, snp_two) in epi_results_dict):
            epi_results_dict[(snp_one, snp_two)] = float(bal_acc)
        
    sorted_epi_results = sorted(epi_results_dict.items(), key=lambda x:x[1], reverse=True)
    #print(sorted_epi_results)

    if(return_dict):
        return sorted_epi_results
    else:
        report_epistasis_results_helper(sorted_epi_results, inter_dis_snps, num_tps, 't-statistic')

def report_mdr_epistasis_results(fname_json:str, inter_dis_snp_indeces:list, fname_epistasis:str,
                                 return_dict:bool = False):
    inter_dis_snps = get_disease_snps_from_epigen_json(fname_json, inter_dis_snp_indeces)
    num_tps = len(inter_dis_snps)
    # print(inter_dis_snps)
    
    epi_results_dict = {}
    sorted_epi_results = {}
    key_lines = []
    with open(fname_epistasis, 'r') as f:
        lines = f.readlines()
        i = 1
        record_lines = False
        while i < len(lines):
            if(lines[i][0:8] == '@Results'):
                if(lines[i+1][0:14] == 'BeginTopModels'):
                    i += 2
                    record_lines = True
            elif(lines[i][0:12] == 'EndTopModels'):
                record_lines = False
                
            if(record_lines):
                key_lines.append(lines[i])

            i += 1
            
    for line in key_lines:
        snp_one = line.split(',')[0]
        part_two = line.split(',')[1]
        snp_two = part_two.split('\t')[0]
        bal_acc = part_two.split('\t')[1]

        if(not (snp_one, snp_two) in epi_results_dict):
            epi_results_dict[(snp_one, snp_two)] = float(bal_acc)
        
    sorted_epi_results = sorted(epi_results_dict.items(), key=lambda x:x[1], reverse=True)
    #print(sorted_epi_results)

    if(return_dict):
        return sorted_epi_results
    else:
        report_epistasis_results_helper(sorted_epi_results, inter_dis_snps, num_tps, 'bal accuracy')

def report_episnp_epistasis_results(fname_json:str, inter_dis_snp_indeces:list, fname_epistasis:str,
                                    return_dict:bool = False, pval_threshold:float = 1e-7):
    inter_dis_snps = get_disease_snps_from_epigen_json(fname_json, inter_dis_snp_indeces)
    num_tps = len(inter_dis_snps)
    # print(inter_dis_snps)
    
    epi_results_dict = {}
    sorted_epi_results = {}
    with open(fname_epistasis, 'r') as f:
        lines = f.readlines()
        i = 0
        while i < len(lines):
            #print(lines[i])
            
            line = lines[i].split(' ')
            line_cleaned = []
            for ele in line:
                if(ele != '' and ele != '\n'):
                    line_cleaned.append(ele)
            
            #print(line_cleaned)

            if(not (line_cleaned[3], line_cleaned[1]) in epi_results_dict):
                metric_value = 0
                if('-' in line_cleaned[6] and ('E-' not in line_cleaned[6])):
                    index = line_cleaned[6].index('-')
                    metric_value = line_cleaned[6][0:index] + 'E-' + line_cleaned[6][(index+1):]
                else:
                    metric_value = line_cleaned[6]
                if(float(metric_value) < pval_threshold):
                    epi_results_dict[(line_cleaned[1], line_cleaned[3])] = float(metric_value)
            i += 1
        sorted_epi_results = sorted(epi_results_dict.items(), key=lambda x:x[1])
        #print(sorted_epi_results)
      
    if(return_dict):
        return sorted_epi_results
    else:
        report_epistasis_results_helper(sorted_epi_results, inter_dis_snps, num_tps, 'p-value')
  
def report_matrix_epistasis_results(fname_json:str, inter_dis_snp_indeces:list, fname_epistasis:str,
                                    return_dict:bool = False):
    
    inter_dis_snps = get_disease_snps_from_epigen_json(fname_json, inter_dis_snp_indeces)
    num_tps = len(inter_dis_snps)
    # print(inter_dis_snps)
    
    epi_results_dict = {}
    sorted_epi_results = {}
    with open(fname_epistasis, 'r') as f:
        lines = f.readlines()
        i = 1
        while i < len(lines):
            line = lines[i].split('\t')

            if(not (line[1], line[0]) in epi_results_dict):
                epi_results_dict[(line[0], line[1])] = float(line[2])
            i += 1
        sorted_epi_results = sorted(epi_results_dict.items(), key=lambda x:x[1])
        #print(sorted_epi_results)
        
    if(return_dict):
        return sorted_epi_results
    else:
        report_epistasis_results_helper(sorted_epi_results, inter_dis_snps, num_tps, 'p-value')

def report_plink_epistasis_results(fname_json:str, inter_dis_snp_indeces:list, fname_epistasis:str,
                                   return_dict = False):
    # Report the following: 
    # Number of true positives identified. 
    # Number of false positives.
    # Average position of true positive(s) when ranked by p-value. 
    
    inter_dis_snps = get_disease_snps_from_epigen_json(fname_json, inter_dis_snp_indeces)
    num_tps = len(inter_dis_snps)
    #print(inter_dis_snps)
    
    # Read in the epistasis epi.qt file
    # Make it a dictionary where in each snp pair maps to a p-value
    # Rank the snp pairs by p-value
    # Identify whether any of the snp pairs are the interacting disease snps (true positive)
    # Count up number of false positives
    # Calculate average position of true positive(s)
    epi_results_dict = {}
    sorted_epi_results = {}
    with open(fname_epistasis, 'r') as f:
        lines = f.readlines()
        i = 1
        while i < len(lines):
            
            line = lines[i].split(' ')
            cleaned_line = []
            for obj in line:
                if(obj != ''):
                    cleaned_line.append(obj)
            #print(cleaned_line[6])
            #print(float(cleaned_line[6]))
            epi_results_dict[(cleaned_line[1], cleaned_line[3])] = float(cleaned_line[6])
            i += 1
        sorted_epi_results = sorted(epi_results_dict.items(), key=lambda x:x[1])
        #print(sorted_epi_results)
        
    if(return_dict):
        return sorted_epi_results
    else:
        report_epistasis_results_helper(sorted_epi_results, inter_dis_snps, num_tps, 'p-value')

def generate_slurm_batch_job_for_plink_epistasis(job_name:str, files:list, epi_thresh:float):
    with open(job_name + '.sh', 'w') as f:
        f.write('#!/bin/sh\n')
        f.write('#SBATCH --job-name=SL_Hail_NA_GWAS\n')
        f.write('#SBATCH --nodes=1\n')
        f.write('#SBATCH --ntasks=1\n')
        f.write('#SBATCH --cpus-per-task=8\n')
        f.write('#SBATCH --mem=32gb\n')
        f.write('#SBATCH --time=240:00:00\n')
        f.write('#SBATCH --partition=shared\n\n')
        f.write('cd $SLURM_SUBMIT_DIR\n')
        f.write('module load plink2/1.90b3.42\n')
        for file in files:
            file = file + '_converted'
            f.write('plink --file ' + file + ' --out ' + file + '\n')
            f.write('plink --bfile ' + file + ' --allow-no-sex --fast-epistasis boost --epi1 ' + str(epi_thresh) + ' --out ' 
                    + file + '_epistasis\n')


def get_disease_snps_from_epigen_json(fname:str, inter_dis_snp_indeces:list):
    disease_snp_names = []
    interacting_disease_snp_names = []
    individual_snps = []
    
    with open(fname, 'r') as f:
        lines = f.readlines()
        
        snps_index = lines[0].find('"snps"')
        disease_snps_index = lines[0].find('"disease_snps"')
        mafs_index = lines[0].find('"mafs"')
        
        # SNP DATA PROCESSING
        snps = lines[0][snps_index + 9:disease_snps_index-3]
        
        temp = ''
        append_flag = False
        #print(genotype)
        for char in snps:
            if(char == '['):
                append_flag = True
                continue
            elif(char == ']'):
                append_flag = False
                individual_snps.append(temp.replace('"', '').split(', '))
                temp = ''
                continue
            if(append_flag):
                temp = temp + char
    
        # REPORT DISEASE SNPS
        disease_snps = lines[0][disease_snps_index+17:mafs_index-3]
        disease_snps = disease_snps.split(', ')
        #print("DISEASE SNPs in file: ", fname)
        #print(disease_snps)
        for ds_snp in disease_snps:
            #print(individual_snps[int(ds_snp)])
            disease_snp_names.append(individual_snps[int(ds_snp)])
            
        for tuple_pair in inter_dis_snp_indeces:
            interacting_disease_snp_names.append((disease_snp_names[tuple_pair[0]][0], disease_snp_names[tuple_pair[1]][0]))
    
    return interacting_disease_snp_names



class TestEpigenCodeBase(unittest.TestCase):  
    # Key Assumptions: Assume that starting EpiGEN json is valid.
    # It contains >= 1 SNPs, >= 1 individuals, valid genotypes, 
    # valid phenotypes, valid snp descriptions, valid number of 
    # disease snps, and valid number of mafs. 
    
    # The tests below heavily rely on following sample json:
    # Impure_Multiplicative_OnePair_BaselineAlpha10_InteractionAlpha2_Chr1_CEU_SNP5_IND10_MAF005_02
    
    def test_extract_key_info_from_EpiGEN_json(self):
        indir = '...'
        fname = indir + 'Impure_Multiplicative_OnePair_BaselineAlpha10_InteractionAlpha2_Chr1_CEU_SNP5_IND10_MAF005_02.json'
        num_inds, num_snps, phenotypes_array, ind_genotypes, ind_snps = extract_key_info_from_EpiGEN_json(fname)
        self.assertEqual(num_inds, 10)
        self.assertEqual(num_snps, 5)
        
        phenotypes_exp = [10.794285182118749, 8.545975042885308, 30.44122919054639, 19.68319380466618,
                          9.690875814819877, 28.19911854711755, 30.504352747819173, 19.036570361202678,
                          10.28237595272995, 10.739893849639815]
        i = 0
        for val in phenotypes_array:
            self.assertEqual(val, phenotypes_exp[i])
            i += 1
            
        ind_genotypes_exp = [['0', '0', '0', '0', '0', '0', '0', '0', '0', '0'], ['0', '0', '1', '0', '0', '1', '1', '0', '0', '0'],
                             ['0', '0', '0', '0', '0', '0', '1', '2', '0', '0'], ['0', '0', '0', '1', '0', '0', '0', '0', '0', '0'],
                             ['0', '0', '0', '0', '0', '0', '0', '0', '0', '0']]
        
        i = 0
        for sublist in ind_genotypes:
            self.assertEqual(sublist, ind_genotypes_exp[i])
            i += 1
        
        ind_snps_exp = [["rs4480304", "chr1", "73472012", "C", "A"],
                        ["rs6700124", "chr1", "92238430", "A", "G"],
                        ["rs4839391", "chr1", "110795513", "G", "A"],
                        ["rs6678558", "chr1", "215349766", "G", "A"],
                        ["rs4268383", "chr1", "246205143", "G", "A"]]
        
        i = 0
        for sublist in ind_snps:
            self.assertEqual(sublist, ind_snps_exp[i])
            i += 1

    def test_extract_key_info_from_multiple_EpiGEN_jsons(self):
        indir = '...'
        fname1 = indir + 'Impure_Multiplicative_OnePair_BaselineAlpha10_InteractionAlpha2_Chr1_CEU_SNP5_IND10_MAF005_02.json'
        fname2 = indir + 'Pure_Multiplicative_OnePair_BaselineAlpha10_InteractionAlpha3_Chr1_CEU_SNP5_IND10_MAF040_045.json'
        fname3 = indir + 'Impure_Multiplicative_OnePair_BaselineAlpha10_InteractionAlpha2_Chr1_CEU_SNP5_IND10_MAF005_02_converted.json'
        
        num_inds, num_snps, phenotypes_array, ind_genotypes, ind_snps = extract_key_info_from_multiple_EpiGEN_jsons([fname1, fname2, fname3])
        num_inds1, num_snps1, phenotypes_array1, ind_genotypes1, ind_snps1 = extract_key_info_from_EpiGEN_json(fname1)
        num_inds2, num_snps2, phenotypes_array2, ind_genotypes2, ind_snps2 = extract_key_info_from_EpiGEN_json(fname2)
        num_inds3, num_snps3, phenotypes_array3, ind_genotypes3, ind_snps3 = extract_key_info_from_EpiGEN_json(fname3)
        
        self.assertEqual(num_inds, num_inds1 + num_inds2 + num_inds3)
        self.assertEqual(num_snps, num_snps1 + num_snps2 + num_snps3)
        
        i = 0
        for val in phenotypes_array1:
            self.assertEqual(val, phenotypes_array[i])
            i += 1
        for val in phenotypes_array2:
            self.assertEqual(val, phenotypes_array[i])
            i += 1
        for val in phenotypes_array3:
            self.assertEqual(val, phenotypes_array[i])
            i += 1
        
        total_genotypes = []
        total_snps = []
        for ind_geno in ind_genotypes1:
            total_genotypes.append(ind_geno)
        for ind_geno in ind_genotypes2:
            total_genotypes.append(ind_geno)
        for ind_geno in ind_genotypes3:
            total_genotypes.append(ind_geno)
            
        i = 0
        for ind_geno in total_genotypes:
            j = 0
            for val in ind_geno:
                self.assertEqual(val, ind_genotypes[i][j])
                j += 1
            i += 1
            
            
        for ind_snp in ind_snps1:
            total_snps.append(ind_snp)
        for ind_snp in ind_snps2:
            total_snps.append(ind_snp)
        for ind_snp in ind_snps3:
            total_snps.append(ind_snp)
            
        i = 0
        for ind_snp in total_snps:
            j = 0
            for val in ind_snp:
                self.assertEqual(val, ind_snps[i][j])
                j += 1
            i += 1

    def test_convert_ped_and_map_to_SNPassoc_csv(self):
        indir = '...'
        fname1 = indir + 'Impure_Multiplicative_OnePair_BaselineAlpha10_InteractionAlpha2_Chr1_CEU_SNP5_IND10_MAF005_02.map'
        fname2 = indir + 'Impure_Multiplicative_OnePair_BaselineAlpha10_InteractionAlpha2_Chr1_CEU_SNP5_IND10_MAF005_02.ped'
        convert_ped_and_map_to_SNPassoc_csv(fname2, fname1)
        
        fname3 = indir + 'Impure_Multiplicative_OnePair_BaselineAlpha10_InteractionAlpha2_Chr1_CEU_SNP5_IND10_MAF005_02_SNPassoc.csv'
        
        snp_names_exp = []
        with open(fname1, 'r') as f:
            lines = f.readlines()
            lines2 = lines[1:]
            for line in lines2:
                snp_names_exp.append(line.split('\t')[1])
                
        
        phenotypes_exp = []
        snp_sets_exp = []
        with open(fname2, 'r') as f:
            lines = f.readlines()
            lines2 = lines[1:]
            for line in lines2:
                temp = line.split('\t')
                phenotype = temp[5]
                snp_values = temp[6:]
                
                phenotypes_exp.append(phenotype)
                snp_sets_exp.append(snp_values)
                
        snp_names = []
        phenotypes = []
        snp_sets = []
        with open(fname3, 'r') as f:
            lines = f.readlines()
            temp = lines[0].split('\t')
            snp_names = temp[2:]
            for line in lines[1:]:
                temp = line.split('\t')
                phenotypes.append(temp[1])
                snp_sets.append(temp[2:])
                
        i = 0
        for snp in snp_names:
            snp = snp.strip()
            self.assertEqual(snp, snp_names_exp[i])
            i += 1
            
        i = 0
        for pheno in phenotypes:
            self.assertEqual(pheno, phenotypes_exp[i])
            i += 1
            
        
        
        snp_sets2 = []
        for snp_set in snp_sets:
            temp = []
            for snp_pair in snp_set:
                snp_pair = snp_pair.strip()
                temp.append(snp_pair[0])
                temp.append(snp_pair[1])
            snp_sets2.append(temp)
        
        i = 0
        for snp_set in snp_sets2:
            j = 0
            for snp in snp_set:
                self.assertEqual(snp, snp_sets_exp[i][j].strip())
                j += 1
                
            i += 1
    
    def test_convert_epigen_json_to_EpiSNP_files(self):
        indir = '...'
        fname = indir + 'Impure_Multiplicative_OnePair_BaselineAlpha10_InteractionAlpha2_Chr1_CEU_SNP5_IND10_MAF005_02.json'
        convert_epigen_json_to_EpiSNP_files(fname)
        
        num_inds, num_snps, phenotypes_array, ind_genotypes, ind_snps = extract_key_info_from_EpiGEN_json(fname)
        
        fname1 = indir + 'Impure_Multiplicative_OnePair_BaselineAlpha10_InteractionAlpha2_Chr1_CEU_SNP5_IND10_MAF005_02.txt'
        fname2 = indir + 'Impure_Multiplicative_OnePair_BaselineAlpha10_InteractionAlpha2_Chr1_CEU_SNP5_IND10_MAF005_02_chr1.dat'
        
        
        snp_names2 = []
        genotypes2 = []
        with open(fname2, 'r') as f:
            lines = f.readlines()
            line = lines[0]
            line = line.split('\t')
            snp_names2 = line[3]
            snp_names2 = snp_names2.split(',')
            
            for line in lines[1:]:
                line = line.split('\t')
                geno = line[3].split(',')
                genotypes2.append(geno)
                
        phenotypes2 = []
        with open(fname1, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                line = line.split('\t')
                phenotypes2.append(line[4])
                
        i = 0
        for snp in ind_snps:
            self.assertEqual(snp[0], snp_names2[i].strip())
            i += 1
            
        i = 0
        for phenotype in phenotypes2:
            self.assertEqual(float(phenotype), phenotypes_array[i])
            i += 1
            
        i = 0
        for ind_geno in ind_genotypes:
            j = 0
            for geno in ind_geno:
                self.assertEqual(geno, genotypes2[j][i].strip())
                j += 1
                
            i += 1
    
    def test_convert_epigen_json_to_genotype_csv(self):
        indir = '...'
        fname = indir + 'Impure_Multiplicative_OnePair_BaselineAlpha10_InteractionAlpha2_Chr1_CEU_SNP5_IND10_MAF005_02.json'
        num_inds, num_snps, phenotypes_array, ind_genotypes, ind_snps = extract_key_info_from_EpiGEN_json(fname)
        convert_epigen_json_to_genotype_csv(fname, True)
        fname1 = indir + 'Impure_Multiplicative_OnePair_BaselineAlpha10_InteractionAlpha2_Chr1_CEU_SNP5_IND10_MAF005_02_QMDR.csv'
        
        snp_names2 = []
        genotypes2 = []
        phenotypes2 = []
        with open(fname1, 'r') as f:
            lines = f.readlines()
            line = lines[0].split(',')
            snp_names2 = line[:-1]
            for line in lines[1:]:
                line = line.split(',')
                genotypes2.append(line[:-1])
                phenotypes2.append(line[-1])
        
        
        i = 0
        for snp in ind_snps:
            self.assertEqual(snp[0], snp_names2[i].strip())
            i += 1
        
        i = 0
        for phenotype in phenotypes2:
            self.assertEqual(float(phenotype), phenotypes_array[i])
            i += 1
        
        i = 0
        for ind_geno in ind_genotypes:
            j = 0
            for geno in ind_geno:
                self.assertEqual(geno, genotypes2[j][i].strip())
                j += 1
                
            i += 1
    
    def test_convert_epigent_qt_to_cc(self):
        indir = '...'
        fname = indir + 'Impure_Multiplicative_OnePair_BaselineAlpha10_InteractionAlpha2_Chr1_CEU_SNP5_IND10_MAF005_02.json'
        convert_epigent_qt_to_cc(fname)
        
        fname2 = indir + 'Impure_Multiplicative_OnePair_BaselineAlpha10_InteractionAlpha2_Chr1_CEU_SNP5_IND10_MAF005_02_converted.json'
        
        num_inds, num_snps, phenotypes_array, ind_genotypes, ind_snps = extract_key_info_from_EpiGEN_json(fname)
        num_inds2, num_snps2, phenotypes_array2, ind_genotypes2, ind_snps2 = extract_key_info_from_EpiGEN_json(fname2)
        
        self.assertEqual(num_inds, num_inds2)
        self.assertEqual(num_snps, num_snps2)
        
        i = 0
        for ind_geno in ind_genotypes:
            j = 0
            for geno in ind_geno:
                self.assertEqual(geno, ind_genotypes2[i][j])
                j += 1
            i += 1
            
        i = 0
        for snp in ind_snps:
            j = 0
            for info in snp:
                self.assertEqual(info, ind_snps2[i][j])
                j += 1
            i += 1
        
        
        # Since this is mult, interaction alpha2 the treshold is 40
        i = 0
        for pheno in phenotypes_array:
            if(pheno > 40):
                temp = 2
            else:
                temp = 1
                
            self.assertEqual(temp, phenotypes_array2[i])
            i += 1
            
        indir = '...'
        fname = indir + 'Impure_Multiplicative_OnePair_BaselineAlpha10_InteractionAlpha2_Chr1_CEU_SNP5_IND10_MAF005_02_MODIFIED.json'
        convert_epigent_qt_to_cc(fname)
        
        fname2 = indir + 'Impure_Multiplicative_OnePair_BaselineAlpha10_InteractionAlpha2_Chr1_CEU_SNP5_IND10_MAF005_02_MODIFIED_converted.json'
    
        num_inds, num_snps, phenotypes_array, ind_genotypes, ind_snps = extract_key_info_from_EpiGEN_json(fname)
        num_inds2, num_snps2, phenotypes_array2, ind_genotypes2, ind_snps2 = extract_key_info_from_EpiGEN_json(fname2)
    
        # Since this is mult, interaction alpha2 the treshold is 40
        i = 0
        for pheno in phenotypes_array:
            if(pheno > 40):
                temp = 2
            else:
                temp = 1
                
            self.assertEqual(temp, phenotypes_array2[i])
            i += 1
    
    def test_convert_epigen_to_ped_and_map(self):
        indir = '...'
        fname = indir + 'Impure_Multiplicative_OnePair_BaselineAlpha10_InteractionAlpha2_Chr1_CEU_SNP5_IND10_MAF005_02.json'
        convert_epigen_to_ped_and_map(fname)
        
        num_inds, num_snps, phenotypes_array, ind_genotypes, ind_snps = extract_key_info_from_EpiGEN_json(fname)
        fname = indir + 'Impure_Multiplicative_OnePair_BaselineAlpha10_InteractionAlpha2_Chr1_CEU_SNP5_IND10_MAF005_02.map'
        
        snp_names2 = []
        with open(fname, 'r') as f:
            lines = f.readlines()
            lines2 = lines[1:]
            for line in lines2:
                snp_names2.append(line.split('\t')[1])
                
        #print(snp_names2)
        num_snps2 = len(snp_names2)
        
        fname = indir + 'Impure_Multiplicative_OnePair_BaselineAlpha10_InteractionAlpha2_Chr1_CEU_SNP5_IND10_MAF005_02.ped'
        
        phenotypes2 = []
        snp_sets2 = []
        with open(fname, 'r') as f:
            lines = f.readlines()
            lines2 = lines[1:]
            for line in lines2:
                temp = line.split('\t')
                phenotype = temp[5]
                snp_values = temp[6:]
                
                phenotypes2.append(phenotype)
                snp_sets2.append(snp_values)
                
        #print(phenotypes2)
        #print(snp_sets2)
        
        num_inds2 = len(snp_sets2)
        
        # Need to convert ind_genotypes into letter format using ind_snps to then compare against snp_sets2
        snp_sets_converted = []
        i = 0
        while i < num_inds2:
            j = 0
            temp = []
            for ind in ind_genotypes:
                val = ind[i]
                if(val == '0'):
                    temp.append(ind_snps[j][3])
                    temp.append(ind_snps[j][3])
                elif(val == '1'):
                    temp.append(ind_snps[j][3])
                    temp.append(ind_snps[j][4])
                elif(val == '2'):
                    temp.append(ind_snps[j][4])
                    temp.append(ind_snps[j][4])
                j += 1
            snp_sets_converted.append(temp)
                
            i += 1
                
        #print(snp_sets_converted)
        
        # Now check that num_snps == num_snps2; num_inds == num_inds2; phenotypes_array == phenotypes2
        # snp_sets2 == snp_sets_converted; snp_names2 match ind_snps
        
        #print(snp_names2)
        #print(ind_snps)
        
        self.assertEqual(num_snps, num_snps2)
        self.assertEqual(num_inds, num_inds2)
        
        i = 0
        for phenotype in phenotypes_array:
            self.assertEqual(phenotype, float(phenotypes2[i]))
            i += 1
            
        i = 0
        for snp_set in snp_sets2:
            j = 0
            for snp in snp_set:
                self.assertEqual(snp[0], snp_sets_converted[i][j])
                j += 1
            i += 1
            
        i = 0
        for snp in ind_snps:
            self.assertEqual(snp[0], snp_names2[i])
            i += 1
    
    
    def test_report_epistasis_results_helper(self):
        sorted_results = [(('rs111111', 'rs222222'), 5e-12), (('rs33333333', 'rs444444444'), 3e-10)]
        inter_dis_snps = [('rs111111', 'rs222222')]
        num_tps = 1
        metric = 'p-value'
        
        one, two, three, four = report_epistasis_results_helper(sorted_results, inter_dis_snps, num_tps, metric, False)
        
        self.assertEqual(one, 1)
        self.assertEqual(two, 1)
        self.assertEqual(three, 1)
        self.assertEqual(four, 1)
        
        
        sorted_results = [(('rs3553635', 'rs45345345'), 5e-25), (('rs543534537676', 'rs788657567567'), 5e-20), 
                          (('rs111111', 'rs222222'), 5e-12), (('rs33333333', 'rs444444444'), 3e-10)]
        inter_dis_snps = [('rs111111', 'rs222222')]
        num_tps = 1
        metric = 'p-value'
        
        one, two, three, four = report_epistasis_results_helper(sorted_results, inter_dis_snps, num_tps, metric, False)
        
        self.assertEqual(one, 1)
        self.assertEqual(two, 3)
        self.assertEqual(three, 3)
        self.assertEqual(four, 3)
        
        sorted_results = [(('rs3553635', 'rs45345345'), 5e-25), (('rs543534537676', 'rs788657567567'), 5e-20), 
                          (('rs111111', 'rs222222'), 5e-12), (('rs33333333', 'rs444444444'), 3e-10)]
        inter_dis_snps = [('rs111111', 'rs222222'), ('rs2423424', 'rs786735345')]
        num_tps = 2
        metric = 'p-value'
        
        one, two, three, four = report_epistasis_results_helper(sorted_results, inter_dis_snps, num_tps, metric, False)
        
        self.assertEqual(one, 1)
        self.assertEqual(two, 3)
        self.assertEqual(three, 3)
        self.assertEqual(four, 4)
    
    def test_report_MIDESP_epistasis_results(self):
        indir = '...'
        file = indir + 'Impure_Dominant_TwoPairs_BaselineAlpha10_InteractionAlpha16_Chr1_CEU_SNP1000_IND1000_MAF005_02.json'
        file2 = indir + 'Impure_Dominant_TwoPairs_BaselineAlpha10_InteractionAlpha16_Chr1_CEU_SNP1000_IND1000_MAF005_02.tped.epiNoAPC'
        
        sorted_results = report_MIDESP_epistasis_results(file, [(0,1),(2,3)], file2, True)
        exp_results = [(('rs12088027', 'rs11163017'), 0.2991390655966821),
                       (('rs12088027', 'rs17021870'), 0.29906200705010216),
                       (('rs12088027', 'rs7518759'), 0.29884142592538815)]
        
        i = 0
        for result in sorted_results:
            self.assertEqual(result, exp_results[i])
            i += 1
        
    def test_report_boost_epistasis_results(self):
        indir = '...'
        file = indir + 'Impure_Multiplicative_OnePair_BaselineAlpha10_InteractionAlpha125_Chr1_CEU_SNP1000_IND1000_MAF005_02_converted.json'
        file2 = indir + 'Impure_Multiplicative_OnePair_BaselineAlpha10_InteractionAlpha125_Chr1_CEU_SNP1000_IND1000_MAF005_02_converted_epistasis.epi.cc'
        
        sorted_results = report_boost_epistasis_results(file, [(0,1)], file2, True)
        exp_results = [(('rs945969', 'rs17462357'), 6.319e-12),
                       (('rs945969', 'rs17119445'), 0.9999),
                       (('rs17119445', 'rs12097141'), 1)]
        
        i = 0
        for result in sorted_results:
            self.assertEqual(result, exp_results[i])
            i += 1
        
    def test_report_qmdr_epistasis_results(self):
        indir = '...'
        file = indir + 'Impure_Dominant_OnePair_BaselineAlpha10_InteractionAlpha8_Chr1_CEU_SNP1000_IND1000_MAF005_02.json'
        file2 = indir + 'Impure_Dominant_OnePair_BaselineAlpha10_InteractionAlpha8_Chr1_CEU_SNP1000_IND1000_MAF005_02_QMDR.csv_analysis.txt'
        
        sorted_results = report_qmdr_epistasis_results(file, [(0,1)], file2, True)
        exp_results = [(('rs12757683', 'rs12565324'), 182.55568),
                       (('rs2710888', 'rs7553373'), 54.254196)]
        
        i = 0
        for result in sorted_results:
            self.assertEqual(result, exp_results[i])
            i += 1
        
    def test_report_mdr_epistasis_results(self):
        indir = '...'
        file = indir + 'Impure_Dominant_OnePair_BaselineAlpha10_InteractionAlpha8_Chr1_CEU_SNP1000_IND1000_MAF005_02.json'
        file2 = indir + 'Impure_Dominant_OnePair_BaselineAlpha10_InteractionAlpha8_Chr1_CEU_SNP1000_IND1000_MAF005_02_converted_QMDR.csv_analysis.txt'
        
        sorted_results = report_mdr_epistasis_results(file, [(0,1)], file2, True)
        exp_results = [(('rs7553373', 'rs6697433'), 0.9403193),
                       (('rs17101563', 'rs7553373'), 0.9370334),
                       (('rs7553373', 'rs4926450'), 0.93621504),
                       (('rs6684002','rs7553373'), 0.9356674),
                       (('rs7553373', 'rs7522427'), 0.93511975)]
        
        i = 0
        for result in sorted_results:
            self.assertEqual(result, exp_results[i])
            i += 1
        
    def test_report_episnp_epistasis_results(self):
        indir = '...'
        file = indir + 'Impure_Dominant_OnePair_BaselineAlpha10_InteractionAlpha8_Chr1_CEU_SNP1000_IND1000_MAF005_02.json'
        file2 = indir + 'Impure_Dominant_OnePair_BaselineAlpha10_InteractionAlpha8_Chr1_CEU_SNP1000_IND1000_MAF005_02_pairwise_net.out'
        
        sorted_results = report_episnp_epistasis_results(file, [(0,1)], file2, True)
        exp_results = [(('rs1002706', 'rs6697433'), 0.128E-16),
                       (('rs873458', 'rs12139261'), 0.705E-11),
                       (('rs10803420', 'rs12139261'), 0.106E-10)]
        
        i = 0
        for result in sorted_results:
            self.assertEqual(result, exp_results[i])
            i += 1
        
    def test_report_matrix_epistasis_results(self):
        pass
        
    def test_report_plink_epistasis_results(self):
        pass
        
    def test_get_disease_snps_from_epigen_json(self):
        indir = '...'
        fname = indir + 'Impure_Multiplicative_OnePair_BaselineAlpha10_InteractionAlpha2_Chr1_CEU_SNP5_IND10_MAF005_02.json'
        inter_dis_snps = get_disease_snps_from_epigen_json(fname, [(0,1)])
        num_inds, num_snps, phenotypes_array, ind_genotypes, ind_snps = extract_key_info_from_EpiGEN_json(fname)
        
        self.assertEqual(inter_dis_snps, [(ind_snps[3][0],ind_snps[0][0])])
        
        indir = '...'
        fname = indir + 'Pure_Multiplicative_OnePair_BaselineAlpha10_InteractionAlpha3_Chr1_CEU_SNP5_IND10_MAF040_045.json'
        inter_dis_snps = get_disease_snps_from_epigen_json(fname, [(0,1)])
        num_inds, num_snps, phenotypes_array, ind_genotypes, ind_snps = extract_key_info_from_EpiGEN_json(fname)
        
        self.assertEqual(inter_dis_snps, [(ind_snps[4][0],ind_snps[2][0])])

def test_suite():
    # Unit Tests
    unit_test_suite = unittest.TestSuite()
    unit_test_suite.addTest(TestEpigenCodeBase('test_extract_key_info_from_EpiGEN_json'))
    unit_test_suite.addTest(TestEpigenCodeBase('test_extract_key_info_from_multiple_EpiGEN_jsons'))
    unit_test_suite.addTest(TestEpigenCodeBase('test_convert_ped_and_map_to_SNPassoc_csv'))
    unit_test_suite.addTest(TestEpigenCodeBase('test_convert_epigen_json_to_EpiSNP_files'))
    unit_test_suite.addTest(TestEpigenCodeBase('test_convert_epigen_json_to_genotype_csv'))
    unit_test_suite.addTest(TestEpigenCodeBase('test_convert_epigent_qt_to_cc'))
    unit_test_suite.addTest(TestEpigenCodeBase('test_convert_epigen_to_ped_and_map'))
    unit_test_suite.addTest(TestEpigenCodeBase('test_report_epistasis_results_helper'))
    unit_test_suite.addTest(TestEpigenCodeBase('test_report_MIDESP_epistasis_results'))
    unit_test_suite.addTest(TestEpigenCodeBase('test_report_boost_epistasis_results'))
    unit_test_suite.addTest(TestEpigenCodeBase('test_report_qmdr_epistasis_results'))
    unit_test_suite.addTest(TestEpigenCodeBase('test_report_mdr_epistasis_results'))
    unit_test_suite.addTest(TestEpigenCodeBase('test_report_episnp_epistasis_results'))
    # unit_test_suite.addTest(TestEpigenCodeBase('test_report_matrix_epistasis_results'))
    # unit_test_suite.addTest(TestEpigenCodeBase('test_report_plink_epistasis_results'))
    unit_test_suite.addTest(TestEpigenCodeBase('test_get_disease_snps_from_epigen_json'))
        
    runner = unittest.TextTestRunner()
    runner.run(unit_test_suite)

# Run the unit tests.
test_suite()

# ANALYZE EPISTASIS RESULTS FROM ALL THE TOOLS SELECTED FOR THE STUDY

for file in (files + files_dominant + files_recessive + files_bonus):
#for file in ['Impure_Dominant_OnePair_BaselineAlpha10_InteractionAlpha8_Chr1_CEU_SNP1000_IND1000_MAF005_02']:
    indir = '...'
    fname = indir + file + '.json' 
    # fname = file + '.fam'
    
    # indir5 = '...'
    # fname5 = indir5 + file + '_MatrixEpistasis.txt'
    
    # indir6 = '...'
    # fname6 = indir6 + file + '_converted_epistasis.epi.cc'
    
    # indir7 = '...'
    # fname7 = indir7 + file + '_pairwise_net.out'
    
    # indir8 = '...'
    # fname8 = indir8 + file + '_converted_QMDR.csv_analysis.txt'
    
    indir9 = '...'
    fname9 = indir9 + file + '_QMDR.csv_analysis.txt'
    
    # indir2 = '...'
    # fname2 = indir2 + file + '_epistasis.epi.qt'
    
    # indir3 = '...'
    # fname3 = indir3 + file + '.tped.epiNoAPC'
    
    # fname4 = indir + file
    
    # REMMA RESULT REPORTING VARIABLES BEGIN
    '''
    fname_xml = '...'
    temp_index = file.find('_Chr1') 
    fname_xml += file[0:temp_index] + '.xml'
    
    remma_sim_name = ''
    if('Impure' in file):
        remma_sim_name += 'impure'
    else:
        remma_sim_name += 'pure'
        
    if('Dominant' in file):
        remma_sim_name += '_dom'
    elif('Recessive' in file):
        remma_sim_name += '_rec'
    elif('Multiplicative' in file):
        remma_sim_name += '_mult'
    elif('XOR' in file):
        remma_sim_name += '_xor'
        
    if('OnePair' in file):
        remma_sim_name += '_onepair'
    elif('TwoPairs' in file):
        remma_sim_name += '_twopairs'
    elif('EightPairs' in file):
        remma_sim_name += '_eightpairs'
        
    if('InteractionAlpha125' in file):
        remma_sim_name += '_interalpha125'
    elif('InteractionAlpha15' in file):
        remma_sim_name += '_interalpha15'
    elif('InteractionAlpha2' in file):
        remma_sim_name += '_interalpha2'
    elif('InteractionAlpha3' in file):
        remma_sim_name += '_interalpha3'
    elif('InteractionAlpha8' in file):
        remma_sim_name += '_interalpha8'
    elif('InteractionAlpha16' in file):
        remma_sim_name += '_interalpha16'
    '''  
    
    # REMMA RESULT REPORTING VARIABLES END
    
    # This is based on EpiGEN model files (.xml) wherein I defined the interacting disease snps. 
    # For now we assume that all interactions are second order only and there is only one, two, or eight pairs.
    if(fname.find('OnePair') != -1):
        inter_dis_snp_indeces = [(0,1)]
    elif(fname.find('TwoPairs') != -1):
        inter_dis_snp_indeces = [(0,1),(2,3)]
    elif(fname.find('EightPairs') != -1):
        inter_dis_snp_indeces = [(0,1),(2,3),(4,5),(6,7),(8,9),(10,11),(12,13),(14,15)]
    
    # if(fname.find('Recessive') != -1):
    #     check_epigen_json_for_detectability(fname, inter_dis_snp_indeces, 'Recessive')
    #     print()
    # elif(fname.find('Dominant') != -1):
    #     check_epigen_json_for_detectability(fname, inter_dis_snp_indeces, 'Dominant')
    #     print()
    # elif(fname.find('Multiplicative') != -1):
    #     check_epigen_json_for_detectability(fname, inter_dis_snp_indeces, 'Multiplicative')
    #     print()
    # elif(fname.find("XOR") != -1):
    #     check_epigen_json_for_detectability(fname, inter_dis_snp_indeces, 'XOR')
    #     print()
    
    #convert_epigent_qt_to_cc(fname)
    #convert_epigen_to_ped_and_map(fname)
    #convert_epigen_json_to_genotype_csv(fname, True)
    #convert_epigen_json_to_EpiSNP_files(fname)
    #print(get_disease_snps_from_epigen_json(fname, inter_dis_snp_indeces))
    #print(file[:-35])
    
    # print(file[:-35])
    # report_boost_epistasis_results(fname, inter_dis_snp_indeces, fname6)
    # report_plink_epistasis_results(fname, inter_dis_snp_indeces, fname2)
    # report_matrix_epistasis_results(fname, inter_dis_snp_indeces, fname5)
    # report_episnp_epistasis_results(fname, inter_dis_snp_indeces, fname7)
    # report_mdr_epistasis_results(fname, inter_dis_snp_indeces, fname8)
    # report_qmdr_epistasis_results(fname, inter_dis_snp_indeces, fname9)
    # report_MIDESP_epistasis_results(fname, inter_dis_snp_indeces, fname3)
    # report_remma_epistasis_results(fname, fname_xml, remma_sim_name)
    # print()
    
    #convert_ped_and_map_to_SNPassoc_csv(fname4 + '.ped', fname4 + '.map')
    
    #convert_epigen_to_DMM(fname, file)
    
    # inter_dis_snps = get_disease_snps_from_epigen_json(fname, inter_dis_snp_indeces)
    # print(inter_dis_snps)
    
    # convert_fam_to_pheno(indir, fname)


# directory = '...'
# raw_file = directory + 'ABCD_202209.updated.nodups.curated.cleaned_indivs_chr11_imputed.raw'
# tfam_file = directory + 'ABCD_202209.updated.nodups.curated.cleaned_indivs_chr11.tfam'
# convert_raw_and_tfam_to_csv_and_pheno(raw_file, tfam_file)

# files_all = files + files_dominant + files_recessive + files_bonus
# generate_slurm_batch_job_for_plink_epistasis("run_plink_epistasis_1000snp_cc_fixedpval", files_all, 1e-7)

# MAKING PLOTS FOR EPISTASIS RESULTS

plt.tight_layout()
fig, ax = plt.subplots()
# BAR CHART1
ax.set_xticks([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18, 19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,
                39,40,41,42,43,44,45])
ax.tick_params(axis='both', labelsize=12)
ax.set_xticklabels(['','MDR','','', '','', '','MIDESP','','', '','', '','QMDR','','', '','', '','Matrix','','',
                    '','', '','Plink','','', '','', '','REMMA','','', '','',  '','BOOST','','', '','', 
                    '','EpiSNP','',''], rotation=45)
# Customize tick parameters
plt.tick_params(axis='x', colors='black', length=0, width=2)
ax.set_ylabel('True Positive %', fontsize = 12)
#print(inspect.signature(ax.set_yticks))
h = plt.bar([0,1,2,3,6,7,8,9,12,13,14,15,18,19,20,21,24,25,26,27,30,31,32,33,36,37,38,39,42,43,44,45],
             [50,50,50,59, 33,66,27,50, 33,50,9.1,9.1, 0,16,18,23, 0,16,18,23, 0,0,23,23, 16,0,0,0, 0,0,0,0], 
             color = ['cornflowerblue', 'royalblue', 'blue', 'darkviolet'], width=0.9)
plt.ylim(0, 100)
plt.legend(h, ['Interaction Alpha 1.25', 'Interaction Alpha 1.5', 'Interaction Alpha 2', 'Interaction Alpha 3'], loc='upper right', fontsize = 12)

# BAR CHART2
#fig.set_size_inches(20, 6, forward=True)
# ax.set_xticks([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22])
# ax.tick_params(axis='both', labelsize=12)
# ax.set_xticklabels(['MDR','', '', 'MIDESP','', '', 'Matrix','', '', 'Plink','', '', 'REMMA','', '', 
#                     'QMDR','', '', 'EpiSNP','', '', 'BOOST',''], rotation=45)
# plt.tick_params(axis='x', colors='black', length=0, width=2)
# ax.set_ylabel('True Positive %', fontsize = 12)
# #print(inspect.signature(ax.set_yticks))
# h = plt.bar([0,1, 3,4, 6,7, 9,10, 12,13, 15,16, 18,19, 21,22],
#              [78.6,41.1, 57.1,19.6, 32.1,17.8, 32.1,17.8, 32.1,17.8, 23.2,5.3, 10.7,3.6, 5.3,12.5], 
#              color = ['red', 'blue'],
#              width = 0.8)
# plt.ylim(0, 100)
# plt.legend(h, ['Pure', 'Impure'], 
#            fontsize = 12)


# #VIOLIN PLOT
# ax.set_yticks([1,2,3])
# ax.set_yticklabels(['MDR XOR', 'MIDESP XOR', 'BOOST XOR'])

# #['QMDR Recessive', 'MDR Recessive', 'REMMA Recessive', 'Plink Recessive', 
# # 'Matrix Recessive','MIDESP Recessive', 'BOOST Recessive', 'EpiSNP Recessive']
# ax.set_ylabel('TP Avg. Position')

# ax.violinplot([[63.75, 121.71, 16.3, 85],
#                [42.71, 24.8, 3.5, 4],
#                [1, 1, 1.5, 1.5],
#                 ], positions = [1,2,3], showmeans = True, 
#                 vert = False)

# ax.scatter([63.75, 121.71, 16.3, 85,
#             42.71, 24.8, 3.5, 4,
#             1, 1, 1.5, 1.5,],
#             [1,1,1,1,
#             2,2,2,2,
#             3,3,3,3])

plt.savefig('filename.png', dpi=600, bbox_inches="tight")




