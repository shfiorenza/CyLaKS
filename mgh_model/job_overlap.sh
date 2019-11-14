#!/bin/bash

# === List of SBATCH arguments ===
#SBATCH --account=ucb-summit-smr	# allocations account
#SBATCH --partition=shas			# Type of node to run on
#SBATCH --qos=normal        		# quality of service/queue
#SBATCH --nodes=2					# number of cluster nodes, abbreviated by -N
#SBATCH --ntasks=48					# number of parallel process
#SBATCH --time=24:00:00				# walltime, abbreviated by -t
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

for C_EFF_TETH in 500 1750 3000 4250
do
	for C_EFF_BIND in 2500 3500 4500 
	do
		for K_OFF_RATIO in 5 10 25 50
		do
			K_OFF_II=$(echo "scale=5; 0.143 / $K_OFF_RATIO" | bc)
			SIM_NAME="slide_scan_"
			SIM_NAME+=$C_EFF_TETH
			SIM_NAME+="_"
			SIM_NAME+=$C_EFF_BIND
			SIM_NAME+="_"
			SIM_NAME+=$K_OFF_RATIO
			TEMP_PARAMS="params_temp_"
			TEMP_PARAMS+=$SIM_NAME
			TEMP_PARAMS+=".yaml"
			cp $PARAM_FILE $TEMP_PARAMS
			yq w -i $TEMP_PARAMS motors.c_eff_tether $C_EFF_TETH
			yq w -i $TEMP_PARAMS xlinks.c_eff_bind $C_EFF_BIND
			yq w -i $TEMP_PARAMS xlinks.k_off_ii $K_OFF_II
			# Run sim for these parameter values
			echo "Running new simulation; name is $SIM_NAME"
			# srun -N 1 -n 1 ./sim $TEMP_PARAMS $SIM_NAME &
		done
	done
done
wait
rm params_temp_slide_scan_*
echo "Scan finished"
