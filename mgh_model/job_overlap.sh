#!/bin/bash

# === List of SBATCH arguments ===
#SBATCH --account=ucb-summit-smr		# allocations account
#SBATCH --partition=smem			# Type of node to run on
#SBATCH --nodes=1				# number of cluster nodes, abbreviated by -N
#SBATCH --ntasks=12 				# number of parallel process
#SBATCH --time=168:00:00			# walltime, abbreviated by -t
#SBATCH --job-name=slideScan		 	# name of job
#SBATCH --output=slideScan_%j.out   		# name of the stdout redirection file
#SBATCH --error=slideScan_%j.err		# name of the stderr redirection file

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

for C_EFF_TETH in {50,500,2500}; do
	for E_SCALE in {0,100,225,300}; do
		E_INT=$(echo "scale=3; $E_SCALE * 0.01" | bc)
		SIM_NAME="coop_slide_scan_"
		SIM_NAME+=$C_EFF_TETH
		SIM_NAME+="_"
		SIM_NAME+=$E_SCALE
		TEMP_PARAMS="params_temp_"
		TEMP_PARAMS+=$C_EFF_TETH
		TEMP_PARAMS+="_"
		TEMP_PARAMS+=$E_SCALE
		TEMP_PARAMS+=".yaml"
		cp $PARAM_FILE $TEMP_PARAMS
		yq w -i $TEMP_PARAMS xlinks.interaction_energy $E_INT
		yq w -i $TEMP_PARAMS motors.c_eff_tether $C_EFF_TETH
		# Run sim for these parameter values
		echo "Running new simulation; name is $SIM_NAME"
		srun -N 1 -n 1 ./sim $TEMP_PARAMS $SIM_NAME &
	done
done
wait
rm params_temp_*
echo "Scan finished"
