#!/bin/bash
echo STARTING PRC1 COOP SCAN
PARAM_FILE="params_endtag.yaml"
echo BASE PARAM FILE is $PARAM_FILE
MT_LENGTH=1000;
MOT_CONC=0.0;
for E_INT in 2.25
do
	for XLINK_CONC in 0.92 1.80 9.20 19.00 28.00 38.00
	do
		FILE_NAME="coop_bindB_"
		FILE_NAME+=$E_INT
		FILE_NAME+="_"
		FILE_NAME+=$XLINK_CONC
		echo RUNNING NEW SIM: filename is $FILE_NAME
		TEMP_PARAMS="params_temp_"
		TEMP_PARAMS+=$FILE_NAME
		TEMP_PARAMS+=".yaml"
		cp $PARAM_FILE $TEMP_PARAMS
		yq w -i $TEMP_PARAMS microtubules.length[0] $MT_LENGTH
		yq w -i $TEMP_PARAMS motors.c_bulk $MOT_CONC
		yq w -i $TEMP_PARAMS xlinks.c_bulk $XLINK_CONC
		yq w -i $TEMP_PARAMS xlinks.interaction_energy $E_INT
		# Run sim for these parameter values
		./sim $TEMP_PARAMS $FILE_NAME &
	done
	wait
done
rm params_temp_coop_bind_*
echo END SCAN
