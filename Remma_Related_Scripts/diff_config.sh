# <description>
# Creates subfolder in config folder and copies relevant files from another config folder's same dataset subfolder. 
#
# Usage:
#  $ bash diff_config.sh [config to copy from] [config to copy to] [dataset identification]

# Create directory in config folder to copy to 
mkdir $2/$2_$3
# Copy all files over
cp $1/$1_$3/* $2/$2_$3
# Remove intermediates and output files
cd $2/$2_$3
rm dom* var* slurm* epi* add* 
# Rename REMMA execution files
mv REMMA_$1_$3\.sh REMMA_$2_$3\.sh
mv REMMA_$1_$3\.py REMMA_$2_$3\.py
