#!/bin/bash
echo STARTING PROCESSIVITY SCAN
PARAM_FILE="params_processivity.yaml"
echo BASE PARAM FILE is $PARAM_FILE

for APP_FORCE in 0 2 4 6 8 10 12
do
	FILE_NAME="processivity_${APP_FORCE}pN"
	echo RUNNING NEW SIM: filename is $FILE_NAME
	TEMP_PARAMS="params_temp_${FILE_NAME}.yaml"
	cp $PARAM_FILE $TEMP_PARAMS
    yq w -i $TEMP_PARAMS motors.applied_force ${APP_FORCE}
	# Run sim for these parameter values
	./sim $TEMP_PARAMS $FILE_NAME &
done
wait
rm params_temp_processivity*
echo END SCAN
