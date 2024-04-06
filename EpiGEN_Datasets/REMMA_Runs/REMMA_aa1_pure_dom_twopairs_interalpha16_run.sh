#!/bin/sh
#SBATCH --job-name=SL_REMMA_aa1_pure_dom_twopairs_interalpha16
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32gb
#SBATCH --time=240:00:00
#SBATCH --partition=shared

cd $SLURM_SUBMIT_DIR
module load use.own
module load python/3.8.3
export PYTHONPATH=/gpfs/home/slistopad/.local/lib/python3.8/site-packages:$PYTHONPATH
python3 REMMA_aa1_pure_dom_twopairs_interalpha16.py
