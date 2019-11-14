#!/bin/bash
echo Starting overlap occupancy scan
PARAM_FILE="params_processivity.yaml"
echo Base parameter file is $PARAM_FILE
MOT_CONC=0.0;
for DIFF_SCALE in 2 2.25 2.5 3
do
	D_II=$(echo "scale=4; 0.131 / $DIFF_SCALE" | bc)
	for I_RUN in 1 2 3
	do
		FILE_NAME="diffu_double"
		FILE_NAME+="_"
		FILE_NAME+=$DIFF_SCALE
		FILE_NAME+="_"
		FILE_NAME+=$I_RUN
		echo RUNNING NEW SIM: filename is $FILE_NAME
		TEMP_PARAMS="params_temp_xlink_diffu"
		TEMP_PARAMS+="_"
		TEMP_PARAMS+=$DIFF_SCALE
		TEMP_PARAMS+="_"
		TEMP_PARAMS+=$I_RUN
		TEMP_PARAMS+=".yaml"
		cp $PARAM_FILE $TEMP_PARAMS
		yq w -i $TEMP_PARAMS motors.c_bulk $MOT_CONC
		yq w -i $TEMP_PARAMS xlinks.diffu_coeff_ii $D_II
		# Run sim for these parameter values
		./sim $TEMP_PARAMS $FILE_NAME &
	done
done
wait
rm params_temp_xlink_diffu*
echo END SCAN
