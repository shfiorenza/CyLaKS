#!/bin/bash
SCAN_NAME="lattice_coop"
echo Starting ${SCAN_NAME} scan
PARAM_FILE="params_processivity.yaml"
echo Base parameter file is ${PARAM_FILE}

MOT_CONC=0.420
XLINK_CONC=0.0
for EMAX_BULK in 1.5 2.0 2.5 
do
	for EMAX_SOLO in 0.75 1.0 1.25
	do
		FILE_NAME="${SCAN_NAME}_${EMAX_BULK}_${EMAX_SOLO}"
		echo Running new sim: file name is ${FILE_NAME}
		TEMP_PARAMS="params_temp_${FILE_NAME}.yaml"
		cp $PARAM_FILE $TEMP_PARAMS
		yq w -i $TEMP_PARAMS motors.c_bulk $MOT_CONC
		yq w -i $TEMP_PARAMS xlinks.c_bulk $XLINK_CONC
		yq w -i $TEMP_PARAMS motors.lattice_coop_Emax_solo $EMAX_SOLO
		yq w -i $TEMP_PARAMS motors.lattice_coop_Emax_bulk $EMAX_BULK
		# Run sim for these parameter values
		./sim $TEMP_PARAMS $FILE_NAME &
	done
done
wait
rm params_temp_${SCAN_NAME}*
echo END SCAN
