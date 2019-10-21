#!/bin/bash

# === List of SBATCH arguments ===
#SBATCH --account=ucb-summit-smr		# allocations account
#SBATCH --partition=shas		    	# Type of node to run on
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

while [ $XLINK_CONC -le $MAX_CONC ]
do
    FILE_NAME="coop_225_"
    FILE_NAME+=$XLINK_CONC
	echo "RUNNING NEW SIM: filename is $FILE_NAME"
	TEMP_PARAMS="params_tempB_"
    TEMP_PARAMS+=$XLINK_CONC
	TEMP_PARAMS+=".yaml"
	cp $PARAM_FILE $TEMP_PARAMS
	yq w -i $TEMP_PARAMS microtubules.length[0] $MT_LENGTH
	yq w -i $TEMP_PARAMS motors.c_bulk $MOT_CONC
	yq w -i $TEMP_PARAMS xlinks.c_bulk $XLINK_CONC
	srun -N 1 -n 1 ./sim $TEMP_PARAMS $FILE_NAME &
#    ./sim $TEMP_PARAMS $FILE_NAME &
    XLINK_CONC=$(($XLINK_CONC + 1))
done
wait
rm params_tempB_*
echo "Scan finished"