#!/bin/bash

# === List of SBATCH arguments ===
#SBATCH --nodes=1			# number of cluster nodes, abbreviated by -N
#SBATCH --time=20:00:00		# walltime, abbreviated by -t
#SBATCH --qos=normal        # quality of service/queue
#SBATCH --job-name=scan		# name of job
#SBATCH --output=scan_%j.out     	# name of the stdout redirection file
#SBATCH --error=scan_%j.err		# name of the stderr redirection file
#SBATCH --ntasks=10		    	# number of parallel process

# === Purge all modules and load needed ones ===
module purge

module load intel
module load gsl
module list
echo "Current variables CC=${CC} and CXX=${CXX}"

MT_LENGTH=1000
STEP_SIZE=200
N_STEPS=$(echo "scale=0; $MT_LENGTH/$STEP_SIZE + 1" | bc)
COORD_SHIFT=0
# === Run the program ===
echo "Running on $(hostname --fqdn)"
for i_mt in {0, 1}; do
	for i_run in {1..$N_STEPS}; do
		PARAM_FILE="params_shift_$(i_mt)_$(COORD_SHIFT).yaml"
		SIM_NAME="shift_$(i_mt)_$(COORD_SHIFT).yaml"
		srun -n 1 ./sim $PARAM_FILE $SIM_NAME &
		COORD_SHIFT+=$(STEP_SIZE)
	done
done

wait
