#!/bin/bash
echo STARTING ENDTAG SCAN
PARAM_FILE="params_processivity.yaml"
echo BASE PARAM FILE is $PARAM_FILE
for CONC_SCALE in {20,50,80,120,220,420}; do
MOT_CONC=$(echo "scale=4; $CONC_SCALE * 0.001" | bc)
	FILE_NAME="processivityG"
	FILE_NAME+="_"
	FILE_NAME+=$CONC_SCALE
	FILE_NAME+="pM"
	echo RUNNING NEW SIM: filename is $FILE_NAME
	TEMP_PARAMS="params_temp_"
	TEMP_PARAMS+=$CONC_SCALE
	TEMP_PARAMS+=".yaml"
	cp $PARAM_FILE $TEMP_PARAMS
	yq w -i $TEMP_PARAMS motors.c_bulk $MOT_CONC
	# Run sim for these parameter values
	./sim $TEMP_PARAMS $FILE_NAME &
done
wait
rm params_temp_*
echo END SCAN
