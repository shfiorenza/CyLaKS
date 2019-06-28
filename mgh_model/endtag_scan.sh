#!/bin/bash
echo STARTING ENDTAG SCAN
PARAM_FILE="params_endtag.yaml"
echo BASE PARAM FILE is $PARAM_FILE
for CONC_SCALE in {2,4}; do
	XLINK_CONC=$(echo "scale=2; $CONC_SCALE * 0.1" | bc)
	for N_SITES in {250,500,750,1000,1250,1750}; do 
		FILE_NAME="EndtagREC"
		FILE_NAME+="_"
		FILE_NAME+=$CONC_SCALE
		FILE_NAME+="_"
		FILE_NAME+=$N_SITES
		echo RUNNING NEW SIM: filename is $FILE_NAME
		TEMP_PARAMS="params_temp_"
		TEMP_PARAMS+=$CONC_SCALE
		TEMP_PARAMS+="_"
		TEMP_PARAMS+=$N_SITES
		TEMP_PARAMS+=".yaml"
		cp $PARAM_FILE $TEMP_PARAMS
		yq w -i $TEMP_PARAMS xlinks.concentration $XLINK_CONC
		yq w -i $TEMP_PARAMS microtubules.length[0] $N_SITES
		# Run sim for these parameter values
		./sim $TEMP_PARAMS $FILE_NAME &
	done
	wait
	rm params_temp_*
done
echo END SCAN
