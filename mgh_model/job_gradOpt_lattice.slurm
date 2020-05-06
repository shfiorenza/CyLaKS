#!/bin/bash

# === List of SBATCH arguments ===
#SBATCH --account=ucb-summit-smr	# allocations account
#SBATCH --partition=shas			# Type of node to run on
#SBATCH --qos=long        		# quality of service/queue
#SBATCH --nodes=1					# number of cluster nodes, abbreviated by -N
#SBATCH --ntasks=6					# number of parallel process
#SBATCH --time=168:00:00				# walltime, abbreviated by -t
#SBATCH --job-name=gradOpt		 	# name of job
#SBATCH --output=job_gradOpt_lattice_%j.out   # name of the stdout redirection file
#SBATCH --error=job_gradOpt_lattice_%j.err	# name of the stderr redirection file

# === Purge all modules and load needed ones ===
module purge
module load gcc/8.2.0
module list
echo "Current variables CC=${CC} and CXX=${CXX}"

# === Run the program ===
echo "Running on $(hostname --fqdn)"
echo "Start of gradient descent optimization scan"

srun python gradOpt_lattice_coop.py
