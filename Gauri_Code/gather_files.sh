# <description>
# Creates directory and gathers REMMA input files in directory. 
#
# Usage:
#  $ bash gather_files.sh [name of data to copy] [name of new directory]

# Create new directory with name of second input 
mkdir $2
# copy bed, fam, bim and pheno files to directory
cp REMMA_Runs/REMMA_DATA/$1.bed $2 
cp REMMA_Runs/REMMA_DATA/$1.fam $2 
cp REMMA_Runs/REMMA_DATA/$1.bim $2 
cp REMMA_Runs/REMMA_DATA/$1_Pheno.ped $2 
