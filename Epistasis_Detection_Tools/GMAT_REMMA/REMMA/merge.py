from glob import glob
from scipy.stats import chi2
import sys

def merge(prefix, p_value = False):
	
	pattern = prefix + '.*[0-9]-*[0-9]'
	
	files = glob(pattern)
	
	if(not files):
		print "No files are grabbed"
		exit()
	else:
		print len(files), " files are grabbed"
	files.sort()
	merge_file = prefix + '.merge'
	fout = open(merge_file, 'w')
	if(p_value):
		fin = open(files[0], 'r')
		line = fin.readline()
		line = line.strip()
		stri = line + ' p_value' + '\n'
		fout.write(stri)
		fin.close()
		for i in files:
			print i
			fin = open(i,'r')
			line = fin.readline()
			line = fin.readline()
			while line:
				line = line.strip()
				arr = line.split()
				p = chi2.sf(float(arr[-1]), 1)
				stri = line + ' ' + str(p) + '\n'
				fout.write(stri)
				line = fin.readline()
			fin.close()
	else:
		fin = open(files[0], 'r')
		line = fin.readline()
		fout.write(line)
		fin.close()
		for i in files:
			print i
			fin = open(i,'r')
			line = fin.readline()
			line = fin.readline()
			while line:
				fout.write(line)
				line = fin.readline()
			fin.close()
	
	fout.close()

if(sys.argv[1] == '--help'):
	print "Quick start: merge --prefix epi --p_val 0/1"
	sys.exit()

print "Read the parameters."
if(len(sys.argv) != 5):
	print "The program need 2 parameters! Please check!"
	sys.exit()

label = 0
for i in range(1, 3):
	if(sys.argv[2*i-1] == '--prefix'):
		res_prefix = sys.argv[2*i]
		label = 1
if(label == 0):
	print "Please provide the parameter --prefix"
	sys.exit()

label = 0
for i in range(1, 3):
	if(sys.argv[2*i-1] == '--p_val'):
		p_val = int(sys.argv[2*i])
		label = 1
if(label == 0):
	print "Please provide the parameter --p_val"
	sys.exit()

if (p_val == 1):
	p_val = True
else:
	p_val = False

merge(res_prefix, p_value = p_val)
