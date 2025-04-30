import sys
from statistics import mean

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
	
		# REPORT DISEASE SNPSGGGGGG
		disease_snps = lines[0][disease_snps_index+17:mafs_index-3]
		disease_snps = disease_snps.split(', ')
		#print("DISEASE SNPs in file: ", fname)GGGGGG
		#print(disease_snps)
		for ds_snp in disease_snps:
			#print(individual_snps[int(ds_snp)])GGGGGG
			disease_snp_names.append(individual_snps[int(ds_snp)])
	
		for tuple_pair in inter_dis_snp_indeces:
			interacting_disease_snp_names.append((disease_snps[tuple_pair[0]], disease_snps[tuple_pair[1]]))
			
		return interacting_disease_snp_names

def report_epistasis_results_helper(sorted_epi_results: list, inter_dis_snps:set, num_tps:int, metric:str,
									verbose:bool = True):
	"""Returns true positive count, false positive count, false positive count, and penalized and unpenalized average true positive position given results of tool sorted by metric, list of interacting snps, and number of true positives.
	
	Code loops through sorted results and checks against list of interacting SNP pairs. 

	Prints results if verbose is True. 

	Code is modified from my mentor's code for reporting results from other tools.
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
		
	if(verbose):
		print(tp_pairs)
		print("Number of true positives: ", true_pos_count)
		print("Number of false positives: ", false_pos_count)
		print("Average true positive position (unpenalized) ranked by " + metric + ": ", mean(tp_indeces))
		print("Average true positive position (penalized) ranked by " + metric + ": ", mean(tp_indeces_penalized))
	
	return true_pos_count, false_pos_count, mean(tp_indeces), mean(tp_indeces_penalized)

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
	
def report_remma_epistasis_results(fname_json:str, fname_xml:str, sim_data:str, return_dict = False, print_results=True):
	""" Given dataset information (json, model, and folder identifier), report results for all five configurations. 

	Return dictionary if return_dict is True, else return statistics. Print results if print_results is True. 
	"""

	# Get positions of interacting SNP pairs
	inter_dis_snp_indeces = extract_interaction_pos(fname_xml)
	# Get SNP numbers from positions (mentor's function)
	inter_dis_snps = get_disease_snps_from_epigen_json(fname_json, inter_dis_snp_indeces)
	# Number of true positives is number of pairs
	num_tps = len(inter_dis_snps)

	config_list = ['aa1', 'aa2', 'aa3', 'ad', 'dd']
	epi_results = {}
	# Iterate through all configurations 
	for config in config_list: 
		# Get path to output file 
		if config == 'aa1' or config == 'aa2' or config == 'aa3':
			path = config+'/'+config+'_'+sim_data+'/epiAA_'+config+'_'+sim_data 
		else:
			path = config+'/'+config+'_'+sim_data+'/epi'+config.upper()+'_'+config+'_'+sim_data 
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
				key = (line[0], line[1])
				# Add pair and p-value to results
				if key in epi_results:
					# Handle overlapping results between configurations by recording lowest p-value
					epi_results[key] = min(epi_results[key], float(line[-1][:-1]))
				else:
					epi_results[key] = float(line[-1][:-1])
	# Sort results by p-value
	sorted_epi_results = sorted(epi_results.items(), key=lambda x:x[1])

	if (return_dict):
		# Return sorted results
		return sorted_epi_results
	else: 
		# Calculate and return statistics 
		a, b, c, d = report_epistasis_results_helper(sorted_epi_results, set(inter_dis_snps), num_tps, 'p-value', print_results)
		return (a, len(inter_dis_snps), b, c, d)


if __name__ == '__main__':
	print(report_remma_epistasis_results('json/Pure_Dominant_OnePair_BaselineAlpha10_InteractionAlpha8_Chr1_CEU_SNP1000_IND1000_MAF005_02.json','EpiGEN_Models/Pure_Dominant_OnePair_BaselineAlpha10_InteractionAlpha8.xml', 'pure_dom_onepair_interalpha8'))	
