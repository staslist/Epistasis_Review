# <description>
# Submits all bash scripts in a directory (including subdirectories) to Slurm. 
#
# Usage:
#  $ bash sbatch_scripts.sh [path to directory]

# Navigate to the base directory
cd $1 || exit

# Loop through subdirectories
for dir in */ ; do
    if [ -d "$dir" ]; then  # Check if directory
			cd "$dir"  # Change to the subdirectory
			# Submit all bash scripts
			for script in *.sh; do
				echo "Submitting script: $script"
				sbatch "$script"
			done

			cd ..  # Go back to the base directory
	fi
done

