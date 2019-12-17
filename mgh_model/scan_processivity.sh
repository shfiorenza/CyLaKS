#!/bin/bash
echo STARTING PROCESSIVITY SCAN
PARAM_FILE="params_processivity.yaml"
echo BASE PARAM FILE is $PARAM_FILE

BASE_SEED=198261346419
XLINK_CONC=0.0

for CONC_SCALE in 20 #80 120 220 420 
do
	MOT_CONC=$(echo "scale=4; $CONC_SCALE * 0.001" | bc)
	for SEED_NO in 0 1 2 3 4 5 6 7 8 9
	do
		SEED=$(( $BASE_SEED + $SEED_NO ))
		FILE_NAME="processivity"
		FILE_NAME+="_"
		FILE_NAME+=$CONC_SCALE
		FILE_NAME+="pM_"
		FILE_NAME+=$SEED_NO
		echo RUNNING NEW SIM: filename is $FILE_NAME
		TEMP_PARAMS="params_temp_"
		TEMP_PARAMS+=$FILE_NAME
		TEMP_PARAMS+=".yaml"
		cp $PARAM_FILE $TEMP_PARAMS
		yq w -i $TEMP_PARAMS seed $SEED
		yq w -i $TEMP_PARAMS motors.c_bulk $MOT_CONC
		yq w -i $TEMP_PARAMS xlinks.c_bulk $XLINK_CONC
		# Run sim for these parameter values
		./sim $TEMP_PARAMS $FILE_NAME &
	done
	wait
	rm params_temp_processivity*
done
echo END SCAN
