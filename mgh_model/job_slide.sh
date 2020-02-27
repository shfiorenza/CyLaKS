#!/bin/bash

# === List of SBATCH arguments ===
#SBATCH --account=ucb-summit-smr	# allocations account
#SBATCH --partition=shas			# Type of node to run on
#SBATCH --qos=long        		# quality of service/queue
#SBATCH --nodes=8					# number of cluster nodes, abbreviated by -N
#SBATCH --ntasks=192					# number of parallel process
#SBATCH --time=96:00:00				# walltime, abbreviated by -t
#SBATCH --job-name=slide		 	# name of job
#SBATCH --output=slideScan_%j.out   # name of the stdout redirection file
#SBATCH --error=slideScan_%j.err	# name of the stderr redirection file

# === Purge all modules and load needed ones ===
module purge
module load gcc/8.2.0
module list
echo "Current variables CC=${CC} and CXX=${CXX}"

# === Run the program ===
echo "Running on $(hostname --fqdn)"
echo "Start of sliding scan"

PARAM_FILE="params_overlap.yaml"
echo "Base parameter file is $PARAM_FILE"

BASE_SEED=198261346419

for SEED_NO in 1 2 3 4 5 6 7 8 9 10 11 12  
do

for C_EFF_TETH in 4500
do
	for C_EFF_BIND in 4500 
	do
		for C_XLINK in 0.2 1.0
		do
			CONC_SCALE=$(echo "scale=0; $C_XLINK * 10" | bc)
			for OVERLAP_LENGTH in 75 150 225 300 375 450 525 600
			do
				SIM_NAME="slide_scan_"
#				SIM_NAME+=$C_EFF_TETH
#				SIM_NAME+="_"
#				SIM_NAME+=$C_EFF_BIND
#				SIM_NAME+="_"
				SIM_NAME+=$CONC_SCALE
				SIM_NAME+="_"
				SIM_NAME+=$OVERLAP_LENGTH
				SIM_NAME+="_"
				SIM_NAME+=$SEED_NO
				TEMP_PARAMS="params_temp_"
				TEMP_PARAMS+=$SIM_NAME
				TEMP_PARAMS+=".yaml"
				SEED=$(( $BASE_SEED + $SEED_NO ))
				cp $PARAM_FILE $TEMP_PARAMS
				yq w -i $TEMP_PARAMS seed $SEED
				yq w -i $TEMP_PARAMS microtubules.length[1] $OVERLAP_LENGTH
				yq w -i $TEMP_PARAMS motors.c_eff_tether $C_EFF_TETH
				yq w -i $TEMP_PARAMS xlinks.c_eff_bind $C_EFF_BIND
				yq w -i $TEMP_PARAMS xlinks.c_bulk $C_XLINK
				# Run sim for these parameter values
				echo "Running new simulation; name is $SIM_NAME"
				srun -N 1 -n 1 ./sim $TEMP_PARAMS $SIM_NAME &
			done
		done
	done
done
done
wait
rm params_temp_slide_scan_*
echo "Scan finished"
