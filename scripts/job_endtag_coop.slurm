#!/bin/bash

# === List of SBATCH arguments ===
#SBATCH --account=ucb-summit-smr	# allocations account
#SBATCH --partition=shas			# Type of node to run on
#SBATCH --qos=condo        		# quality of service/queue
#SBATCH --nodes=5					# number of cluster nodes, abbreviated by -N
#SBATCH --ntasks=120					# number of parallel process
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00				# walltime, abbreviated by -t
#SBATCH --job-name=endtags		 	# name of job
#SBATCH --output=run_endtags_%j.out   # name of the stdout redirection file
#SBATCH --error=run_endtags_%j.err	# name of the stderr redirection file

# === Purge all modules and load needed ones ===
module purge
module load gcc/8.2.0
module list

# === Run the program ===
export LD_LIBRARY_PATH=/projects/shfi4480/overlap_analysis/mgh_model/libs/gsl/lib:$LD_LIBRARY_PATH
echo "LD_LIBRARY_PATH set to ${LD_LIBRARY_PATH}"
echo "Current variables CC=${CC} and CXX=${CXX}"
echo "Running on $(hostname --fqdn)"

SCAN_NAME="endtag"
echo Starting ${SCAN_NAME} scan
PARAM_FILE="params_endtag.yaml"
echo Base parameter file is ${PARAM_FILE}

BASE_SEED=198261346419

for N_SITES in 250 500 750 1000 1250 1750
do
	for RANGE in 10 50 100 500 1000
	do
		for SEED_NO in 0 1 2 3
		do
			FILE_NAME="${SCAN_NAME}_${N_SITES}_${RANGE}_${SEED_NO}"
			TEMP_PARAMS="params_temp_${FILE_NAME}.yaml"
			cp $PARAM_FILE $TEMP_PARAMS
			yq w -i ${TEMP_PARAMS} seed $(( $BASE_SEED + $SEED_NO ))
			yq w -i ${TEMP_PARAMS} motors.gaussian_range ${RANGE}
			yq w -i ${TEMP_PARAMS} filaments.n_sites[0] ${N_SITES}
			# Run sim for these parameter values
			echo Running new sim: file name is ${FILE_NAME}
			./sim $TEMP_PARAMS $FILE_NAME &
		done
	done
done
wait
rm params_temp_${SCAN_NAME}*
echo END SCAN