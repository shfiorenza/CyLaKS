#!/bin/bash
echo Starting overlap occupancy scan
PARAM_FILE="params_overlap.yaml"
echo Base parameter file is $PARAM_FILE
MOT_CONC=0.0;
for K_OFF_SCALE in 2.25 5 10
do
	K_OFF=$(echo "scale=3; 0.143 / $K_OFF_SCALE" | bc)
	for C_EFF in 500 1500 2500 5000
	do
		FILE_NAME="overlap_occu"
		FILE_NAME+="_"
		FILE_NAME+=$K_OFF_SCALE
		FILE_NAME+="_"
		FILE_NAME+=$C_EFF
		echo RUNNING NEW SIM: filename is $FILE_NAME
		TEMP_PARAMS="params_temp_overlap_occu"
		TEMP_PARAMS+="_"
		TEMP_PARAMS+=$K_OFF_SCALE
		TEMP_PARAMS+="_"
		TEMP_PARAMS+=$C_EFF
		TEMP_PARAMS+=".yaml"
		cp $PARAM_FILE $TEMP_PARAMS
		yq w -i $TEMP_PARAMS motors.c_bulk $MOT_CONC
		yq w -i $TEMP_PARAMS xlinks.c_eff_bind $C_EFF
		yq w -i $TEMP_PARAMS xlinks.k_off_ii $K_OFF
		# Run sim for these parameter values
		./sim $TEMP_PARAMS $FILE_NAME &
	done
done
wait
rm params_temp_overlap_occu*
echo END SCAN
