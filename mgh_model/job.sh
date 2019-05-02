#!/bin/bash

# === List of SBATCH arguments ===
#SBATCH --account=ucb-summit-smr	# allocations account
#SBATCH --nodes=1				# number of cluster nodes, abbreviated by -N
#SBATCH --ntasks=7				# number of parallel process
#SBATCH --time=15:00:00			# walltime, abbreviated by -t
#SBATCH --qos=normal        	# quality of service/queue
#SBATCH --job-name=scan			# name of job
#SBATCH --output=scan_%j.out    # name of the stdout redirection file
#SBATCH --error=scan_%j.err		# name of the stderr redirection file

# === Purge all modules and load needed ones ===
module purge
module load intel
module load gsl
module list
echo "Current variables CC=${CC} and CXX=${CXX}"

# === Run the program ===
echo "Running on $(hostname --fqdn)"

N_SEEDS=1;
I_SEED=0;
while [ $I_SEED -lt $N_SEEDS ]; do
	for INITIAL_OVERLAP in {50,100,200,400,600,800,1000}; do
		PARAM_FILE="params_slide_"
		PARAM_FILE+=$I_SEED
		PARAM_FILE+="_"
		PARAM_FILE+=$INITIAL_OVERLAP
		PARAM_FILE+=".yaml"
		echo PARAM FILE is $PARAM_FILE
		SIM_NAME="slide_"
		SIM_NAME+=$I_SEED
		SIM_NAME+="_"
		SIM_NAME+=$INITIAL_OVERLAP
		echo SIM NAME is $SIM_NAME
		srun -N 1 -n 1 ./sim $PARAM_FILE $SIM_NAME &
	done
	let I_SEED+=1;
done
wait
