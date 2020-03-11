#!/bin/bash
echo STARTING PROCESSIVITY SCAN
PARAM_FILE="params_processivity.yaml"
echo BASE PARAM FILE is $PARAM_FILE

XLINK_CONC=0.0
for CONC_SCALE in 20 80 120 220 420 
do
	MOT_CONC=$(echo "scale=4; $CONC_SCALE * 0.001" | bc)
	FILE_NAME="processivity_${CONC_SCALE}pM"
	echo RUNNING NEW SIM: filename is $FILE_NAME
	TEMP_PARAMS="params_temp_${FILE_NAME}.yaml"
	cp $PARAM_FILE $TEMP_PARAMS
	yq w -i $TEMP_PARAMS motors.c_bulk $MOT_CONC
	yq w -i $TEMP_PARAMS xlinks.c_bulk $XLINK_CONC
	# Run sim for these parameter values
	./sim $TEMP_PARAMS $FILE_NAME &
done
wait
rm params_temp_processivity*
echo END SCAN
