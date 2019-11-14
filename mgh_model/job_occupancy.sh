#!/bin/bash

# === List of SBATCH arguments ===
#SBATCH --account=ucb-summit-smr		# allocations account
#SBATCH --partition=shas		    	# Type of node to run on
#SBATCH --qos=normal    		    	# quality of service/queue
#SBATCH --nodes=9			        	# number of cluster nodes, abbreviated by -N
#SBATCH --ntasks=216 				    # number of parallel process
#SBATCH --time=24:00:00		    	    # walltime, abbreviated by -t
#SBATCH --job-name=occu    			 	# name of job
#SBATCH --output=occupancyScan_%j.out   # name of the stdout redirection file
#SBATCH --error=occupancyScan_%j.err	# name of the stderr redirection file

# === Purge all modules and load needed ones ===
module purge
module load gcc/8.2.0
module list
echo "Current variables CC=${CC} and CXX=${CXX}"

# === Run the program ===
echo "Running on $(hostname --fqdn)"
echo "Start of coop binding scan"

PARAM_FILE="params_overlap.yaml"
echo BASE PARAM FILE is $PARAM_FILE
BOT_MT_LENGTH=1500;
TOP_MT_LENGTH=750;
MOT_CONC=0.0;
C_EFF_BIND=1400
MAX_C_EFF=5000;
STEP_SIZE=200;

while [ $C_EFF_BIND -le $MAX_C_EFF ]
do
	for SEED in 1 2 3 4 5 6 7 8 9 10 11 12
	do
		FILE_NAME="occupancy_scan_"
		FILE_NAME+=$C_EFF_BIND
		FILE_NAME+="_"
		FILE_NAME+=$SEED
		echo "RUNNING NEW SIM: filename is $FILE_NAME"
		TEMP_PARAMS="params_temp_"
		TEMP_PARAMS+=$FILE_NAME
		TEMP_PARAMS+=".yaml"
		echo "                 paramfile is $TEMP_PARAMS"
		cp $PARAM_FILE $TEMP_PARAMS
		yq w -i $TEMP_PARAMS microtubules.length[0] $BOT_MT_LENGTH
		yq w -i $TEMP_PARAMS microtubules.length[1] $TOP_MT_LENGTH
		yq w -i $TEMP_PARAMS motors.c_bulk $MOT_CONC
		yq w -i $TEMP_PARAMS xlinks.c_eff_bind $C_EFF_BIND
		srun -N 1 -n 1 ./sim $TEMP_PARAMS $FILE_NAME &
	done
	let "C_EFF_BIND+=STEP_SIZE"
done
wait
rm params_temp_occupancy_scan_*
echo "Scan finished"