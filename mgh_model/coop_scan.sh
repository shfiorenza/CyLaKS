#!/bin/bash
echo STARTING ENDTAG SCAN
PARAM_FILE="params_tuning.yaml"
echo BASE PARAM FILE is $PARAM_FILE
#for E_SCALE in {1,5,10,15,20}; do
E_SCALE=40;
	E_INT=$(echo "scale=3; $E_SCALE * 0.1" | bc)
	for CONC_SCALE in {11,21,106,218,325,431}; do
		XLINK_CONC=$(echo "scale=3; $CONC_SCALE * 0.1" | bc)
		FILE_NAME="coopC"
		FILE_NAME+="_"
		FILE_NAME+=$E_SCALE
		FILE_NAME+="_"
		FILE_NAME+=$CONC_SCALE
		echo RUNNING NEW SIM: filename is $FILE_NAME
		TEMP_PARAMS="params_temp_"
		TEMP_PARAMS+=$CONC_SCALE
		TEMP_PARAMS+=".yaml"
		cp $PARAM_FILE $TEMP_PARAMS
		yq w -i $TEMP_PARAMS xlinks.c_bulk $XLINK_CONC
		yq w -i $TEMP_PARAMS xlinks.interaction_energy $E_INT
		# Run sim for these parameter values
		./sim $TEMP_PARAMS $FILE_NAME &
	done
	wait
	rm params_temp_*
#done
echo END SCAN
