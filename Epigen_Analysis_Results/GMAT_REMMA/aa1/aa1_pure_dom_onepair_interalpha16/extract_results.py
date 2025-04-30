import pandas as pd

def report_remma_epistasis_results(fname_json:str):
	df = pd.read_csv('epiAA_aa1_pure_dom_onepair_interalpha16', sep=' ')
	df.drop(columns=['eff', 'chi'])
	df['p_val'] = pd.to_numeric(df['p_val'], downcast='float')
	print(df.head(10))

if __name__ == '__main__':
	report_remma_epistasis_results()	
