#!/bin/bash

# === List of SBATCH arguments ===
#SBATCH --account=ucb-summit-smr		# allocations account
#SBATCH --partition=shas		    	# Type of node to run on
#SBATCH --qos=normal        	# quality of service/queue
#SBATCH --nodes=2				        # number of cluster nodes, abbreviated by -N
#SBATCH --ntasks=48 				    # number of parallel process
#SBATCH --time=24:00:00		    	    # walltime, abbreviated by -t
#SBATCH --job-name=coop     		 	# name of job
#SBATCH --output=coopScan_%j.out   	# name of the stdout redirection file
#SBATCH --error=coopScan_%j.err		# name of the stderr redirection file

# === Purge all modules and load needed ones ===
module purge
module load gcc/8.2.0
module list
echo "Current variables CC=${CC} and CXX=${CXX}"

# === Run the program ===
echo "Running on $(hostname --fqdn)"
echo "Start of coop binding scan"

PARAM_FILE="params_endtag.yaml"
echo BASE PARAM FILE is $PARAM_FILE
MT_LENGTH=1500;
MOT_CONC=0.0;
XLINK_CONC=1;
MAX_CONC=41;

#while [ $XLINK_CONC -le $MAX_CONC ]
for E_SCALE in 15 30 45 60 75 90 105 120
do
	E_INT=$(echo "scale=3; $E_SCALE * 0.01" | bc)
	for XLINK_CONC in 1 2 9 19 28 38
	do
		FILE_NAME="coop_scan"
		FILE_NAME+="_"
		FILE_NAME+=$E_SCALE
		FILE_NAME+="_"
		FILE_NAME+=$XLINK_CONC
		echo "RUNNING NEW SIM: filename is $FILE_NAME"
		TEMP_PARAMS="params_temp_coop"
		TEMP_PARAMS+="_"
		TEMP_PARAMS+=$E_SCALE
		TEMP_PARAMS+="_"
		TEMP_PARAMS+=$XLINK_CONC
		TEMP_PARAMS+=".yaml"
		echo "                 paramfile is $TEMP_PARAMS"
		cp $PARAM_FILE $TEMP_PARAMS
		yq w -i $TEMP_PARAMS microtubules.length[0] $MT_LENGTH
		yq w -i $TEMP_PARAMS motors.c_bulk $MOT_CONC
		yq w -i $TEMP_PARAMS xlinks.c_bulk $XLINK_CONC
		yq w -i $TEMP_PARAMS xlinks.interaction_energy $E_INT
#		srun -N 1 -n 1 ./sim $TEMP_PARAMS $FILE_NAME &
#    	./sim $TEMP_PARAMS $FILE_NAME &
	done
done
wait
rm params_temp_coop_*
echo "Scan finished"