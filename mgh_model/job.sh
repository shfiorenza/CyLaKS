#!/bin/bash

# === List of SBATCH arguments ===
#SBATCH --nodes=1		# number of cluster nodes, abbreviated by -N
#SBATCH --time=20:00:00		# walltime, abbreviated by -t
#SBATCH --qos=normal        	# quality of service/queue
#SBATCH --job-name=mpi_sim	# name of job
#SBATCH --output=sim.out     	# name of the stdout redirection file
#SBATCH --error=sim.err		# name of the stderr redirection file
#SBATCH --ntasks=24         	# number of parallel process

# === Purge all modules and load needed ones ===
module purge

module load intel
module load impi
module load gsl
module list
echo "Current variables CC=${CC} and CXX=${CXX}"

# === Run the program ===
echo "Running on $(hostname --fqdn)"

export OMP_NUM_THREADS=8
mpirun -n 1 ./sim params.yaml run_8
