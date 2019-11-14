#!/bin/bash
echo STARTING PRC1 COOP SCAN
PARAM_FILE="params_endtag.yaml"
echo BASE PARAM FILE is $PARAM_FILE
MT_LENGTH=1000;
MOT_CONC=0.0;
for E_SCALE in 225
do
	E_INT=$(echo "scale=3; $E_SCALE * 0.01" | bc)
	for CONC_SCALE in 9 18 92 190 280 380
	do
		XLINK_CONC=$(echo "scale=3; $CONC_SCALE * 0.1" | bc)
		FILE_NAME="coop_bind_longB"
		FILE_NAME+="_"
		FILE_NAME+=$E_SCALE
		FILE_NAME+="_"
		FILE_NAME+=$CONC_SCALE
		echo RUNNING NEW SIM: filename is $FILE_NAME
		TEMP_PARAMS="params_temp_"
		TEMP_PARAMS+=$E_SCALE
		TEMP_PARAMS+="_"
		TEMP_PARAMS+=$CONC_SCALE
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
rm params_temp_*
echo END SCAN
