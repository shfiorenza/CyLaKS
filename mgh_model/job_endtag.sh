#!/bin/bash

# === List of SBATCH arguments ===
#SBATCH --account=ucb-summit-smr	# allocations account
#SBATCH --partition=shas			# Type of node to run on
#SBATCH --qos=normal        	# quality of service/queue
#SBATCH --nodes=1				# number of cluster nodes, abbreviated by -N
#SBATCH --ntasks=18				# number of parallel process
#SBATCH --time=24:00:00			# walltime, abbreviated by -t
#SBATCH --job-name=endtag			# name of job
#SBATCH --output=endtagScan_%j.out    # name of the stdout redirection file
#SBATCH --error=endtagScan_%j.err		# name of the stderr redirection file

# === Purge all modules and load needed ones ===
module purge
module load gcc/8.2.0
module list
echo "Current variables CC=${CC} and CXX=${CXX}"

# === Run the program ===
echo "Running on $(hostname --fqdn)"
echo "Start of endtag scan"

PARAM_FILE="params_endtag.yaml"
echo "Base parameter file is $PARAM_FILE"

for SEED in 1 2 3 4
do
	for XLINK_CONC_SCALE in 0 1 4
	do
		XLINK_CONC=$(echo "scale=3; $XLINK_CONC_SCALE * 0.1" | bc)
		for MT_LENGTH in 250 500 750 1000 1250 1750
		do
			SIM_NAME="endtag_scan_"
			SIM_NAME+=$XLINK_CONC_SCALE
			SIM_NAME+="_"
			SIM_NAME+=$MT_LENGTH
			SIM_NAME+="_"
			SIM_NAME+=$SEED
			TEMP_PARAMS="params_temp_"
			TEMP_PARAMS+=$SIM_NAME
			TEMP_PARAMS+=".yaml"
			cp $PARAM_FILE $TEMP_PARAMS
			yq w -i $TEMP_PARAMS seed $SEED
			yq w -i $TEMP_PARAMS xlinks.c_bulk $XLINK_CONC
			yq w -i $TEMP_PARAMS microtubules.length[0] $MT_LENGTH
			echo "Running new simulation; name is $SIM_NAME"
			srun -N 1 -n 1 ./sim $TEMP_PARAMS $SIM_NAME &
		done
	done
done
wait
rm params_temp_endtag_scan_*
