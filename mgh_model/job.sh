#!/bin/bash

# === List of SBATCH arguments ===
#SBATCH --nodes=5			# number of cluster nodes, abbreviated by -N
#SBATCH --time=10:00:00		# walltime, abbreviated by -t
#SBATCH --qos=normal        	# quality of service/queue
#SBATCH --job-name=scan			# name of job
#SBATCH --output=scan_%j.out    # name of the stdout redirection file
#SBATCH --error=scan_%j.err		# name of the stderr redirection file
#SBATCH --ntasks=120		    # number of parallel process

# === Purge all modules and load needed ones ===
module purge

module load intel
module load gsl
module list
echo "Current variables CC=${CC} and CXX=${CXX}"

N_SEEDS=24;
MT_LENGTH=1000;
STEP_SIZE=200;
# === Run the program ===
echo "Running on $(hostname --fqdn)"
I_SEED=0;
while [ $I_SEED -lt $N_SEEDS ]; do
	COORD_SHIFT=0;
	while [ $COORD_SHIFT -lt $MT_LENGTH ]; do
		PARAM_FILE="params_shift_stall_"
		PARAM_FILE+=$I_SEED
		PARAM_FILE+="_"
		PARAM_FILE+=$COORD_SHIFT
		PARAM_FILE+=".yaml"
		echo PARAM FILE is $PARAM_FILE
		SIM_NAME="shift_stall_"
		SIM_NAME+=$I_SEED
		SIM_NAME+="_"
		SIM_NAME+=$COORD_SHIFT
		echo SIM NAME is $SIM_NAME
		srun -N 1 -n 1 ./sim $PARAM_FILE $SIM_NAME &
		let COORD_SHIFT+=$STEP_SIZE;
	done
	let I_SEED+=1;
done
wait
