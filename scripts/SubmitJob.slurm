#!/bin/bash

# === List of SBATCH arguments ===
#SBATCH --account=ucb-summit-smr	# allocations account
#SBATCH --partition=shas			# Type of node to run on
#SBATCH --qos=long        		    # quality of service/queue
#SBATCH --nodes=1					# number of cluster nodes, abbreviated by -N
#SBATCH --ntasks=24					# number of parallel process
#SBATCH --cpus-per-task=1           # cpu allocation per process
#SBATCH --time=168:00:00			# walltime, abbreviated by -t
#SBATCH --job-name=shep		 	    # name of job
#SBATCH --output=run_shep_%j.out    # name of the stdout redirection file
#SBATCH --error=run_shep_%j.err	    # name of the stderr redirection file

# === Purge all modules and load needed ones ===
module purge
module load singularity

# === Run the program ===
echo "Current variables CC=${CC} and CXX=${CXX}"
echo "Running on $(hostname --fqdn)"

bash scan_shepherding.sh