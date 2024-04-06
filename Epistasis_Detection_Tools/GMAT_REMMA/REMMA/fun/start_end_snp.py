import numpy as np

def start_end_snp(num_SNP, total_part, which_part):
	temp_test = 0
	temp_which_part = 0
	start_SNP = 0
	end_SNP = 0
	total_test = num_SNP * (num_SNP - 1)/2
	each_test = np.int(total_test/total_part)
	have_SNP = 0
	for i in range(num_SNP):
		temp_test += num_SNP - i - 1
		if temp_test >= each_test:
			start_SNP = end_SNP
			end_SNP = i + 1
			temp_test = 0
			temp_which_part += 1
		
		if temp_which_part == which_part:
			have_SNP = 1
			break
	
	if have_SNP == 0:
		start_SNP = end_SNP
		end_SNP = num_SNP - 1
	
	return start_SNP, end_SNP

