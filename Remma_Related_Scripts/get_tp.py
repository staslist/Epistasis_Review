import extract_results
import statistics
from statistics import mean, median
import os 

def parse_files():
	""" Generates dictionaries containing relevant files for each type of interaction. """

	# Create dictionary of models
	models = {"dom": [], "mult": [], "rec": [], "xor": []}
	
	# Iterate through all models
	for f in os.listdir("./EpiGEN_Models"):
		# Add file to relevant interaction type's list
		if "Dominant" in f:
			models["dom"].append(f)
		elif "Multiplicative" in f:
			models["mult"].append(f)
		elif "Recessive" in f:
			models["rec"].append(f)
		elif "XOR" in f:
			models["xor"].append(f)
	
	# Sort lists in dictionary 
	for key in models.keys():
		models[key] = sorted(models[key])
	
	# Repeat process for all json files
	jsons = {"dom": [], "mult": [], "rec": [], "xor": []}
	
	for f in os.listdir("./json"):
		if "Dominant" in f:
			jsons["dom"].append(f)
		elif "Multiplicative" in f:
			jsons["mult"].append(f)
		elif "Recessive" in f:
			jsons["rec"].append(f)
		elif "XOR" in f:
			jsons["xor"].append(f)
	
	for key in jsons.keys():
		jsons[key] = sorted(jsons[key])
	
	# Repeat process for all dataset identifiers 
	# Dataset identifiers have the format [purity]_[interaction type]_[number of pairs]_[interaction alpha]
	# ex. pure_dom_onepair_interalpha8
	simdata = {"dom": [], "mult": [], "rec": [], "xor": []}
	
	for d in os.listdir("./aa1"):
		if "dom" in d:
			simdata["dom"].append("_".join(d.split("_")[1:]))
		elif "mult" in d:
			simdata["mult"].append("_".join(d.split("_")[1:]))
		elif "rec" in d:
			simdata["rec"].append("_".join(d.split("_")[1:]))
		elif "xor" in d:
			simdata["xor"].append("_".join(d.split("_")[1:]))
	
	for key in simdata.keys():
		simdata[key] = sorted(simdata[key])
	
	return jsons, models, simdata

def epi_type_analysis():
	""" Print all REMMA results statistics based on interaction type. """

	jsons, models, simdata = parse_files()
	print("*EPISTASIS TYPE ANALYSIS*")
	# Iterate through interacion types
	for key in models.keys():
		total_tp = 0;
		total_p = 0
		total_fp = 0;
		unpen_tp_pos = []
		pen_tp_pos = []
		
		# Get results for all datasets in type
		for i in range(len(models[key])):
			results = extract_results.report_remma_epistasis_results('json/'+jsons[key][i], 'EpiGEN_Models/'+models[key][i], simdata[key][i], print_results=False)
			total_tp += results[0]
			total_p += results[1]
			total_fp += results[2]
			unpen_tp_pos.append(results[3])
			pen_tp_pos.append(results[4])

		# Print results 
		print(key+"--------------------------------------------------------------") 
		print("True positives: "+str(total_tp)+"/"+str(total_p))
		if total_tp == 0:
			print()
			continue
		unpen_tp_pos_drop_zeros = [x for x in unpen_tp_pos if x != 0] 
		pen_tp_pos_drop_zeros = [x for x in pen_tp_pos if x != 0] 
		print("Unpenalized TP positions (no zeros): "+str(unpen_tp_pos_drop_zeros))
		if len(unpen_tp_pos_drop_zeros) != 0:
			print("Unpenalized mean: "+str(mean(unpen_tp_pos_drop_zeros)))
			print("Unpenalized median: "+str(median(unpen_tp_pos_drop_zeros)))
		else:
			print("Unpenalized mean and median: NA")
		print("Penalized TP positions (no zeros): "+str(pen_tp_pos_drop_zeros))
		if len(pen_tp_pos_drop_zeros) != 0:
			print("Penalized mean: "+str(mean(pen_tp_pos_drop_zeros)))
		else:
			print("Unpenalized mean and median: NA")
		print()

def interalpha_purity_analysis():
	""" Print all REMMA results statistics based on purity and interaction alpha. """


	jsons, models, simdata = parse_files()

	# Generate new dictionary with interaction alphas and purity for each interaction type
	epi_type = {"mult":{"125": [], "15": [], "2": [], "3": [], "Pure":[], "Impure":[]}, "dom":{"8":[], "16":[], "Pure":[], "Impure":[]}, "rec":{"8":[], "16":[], "Pure":[], "Impure":[]}, "xor":{"8":[], "16":[], "Pure":[], "Impure":[]}} 
	
	for t in epi_type.keys():
		for i in range(len(simdata[t])):
			if simdata[t][i][0] == "i":
				epi_type[t]["Impure"].append((jsons[t][i], models[t][i], simdata[t][i]))
			else:
				epi_type[t]["Pure"].append((jsons[t][i], models[t][i], simdata[t][i]))
			epi_type[t][simdata[t][i].split("a")[-1]].append((jsons[t][i], models[t][i], simdata[t][i])) 
	
	print("*INTERALPHA/PURITY RESULTS*")

	# Print relevant results for each purity and interaction alpha type 
	for key in epi_type.keys():
		print(key+"--------------------------------------------------------------")
		for setting in epi_type[key].keys():
			print(setting+":")
			total_tp = 0;
			total_p = 0
			total_fp = 0;
			pen_tp_pos = []
			
			for case in epi_type[key][setting]:
				results = extract_results.report_remma_epistasis_results('json/'+case[0], 'EpiGEN_Models/'+case[1], case[2], print_results=False)
				total_tp += results[0]
				total_p += results[1]
				total_fp += results[2]
				pen_tp_pos.append(results[4])
			print("True positives: "+str(total_tp)+"/"+str(total_p))
			print("Penalized TP positions (no zeros): "+str(pen_tp_pos))
			print()

if __name__=="__main__":
	#epi_type_analysis()
	interalpha_purity_analysis()
