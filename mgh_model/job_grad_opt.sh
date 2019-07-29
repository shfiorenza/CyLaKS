#!/bin/bash

# === List of SBATCH arguments ===
#SBATCH --account=ucb-summit-smr		# allocations account
#SBATCH --partition=smem		# node type
#SBATCH --qos=condo        		# quality of service/queue
#SBATCH --nodes=1				# number of cluster nodes, abbrv. by -N
#SBATCH --ntasks=12				# number of parallel process
#SBATCH --time=7-00:00:00		# walltime, abbreviated by -t
#SBATCH --job-name=least_squares		# name of job
#SBATCH --output=least_squares-%j.out   # name of the stdout redirection file
#SBATCH --error=least_squares-%j.err	# name of the stderr redirection file

# === Purge all modules and load needed ones ===
module purge
module load gcc
module load gsl
module load python
module load matlab
module list
echo "Current variables CC=${CC} and CXX=${CXX}"

# === Run the program ===
echo "Running on $(hostname --fqdn)"

srun -N 1 -n 12 ./gradient_optimization.py
