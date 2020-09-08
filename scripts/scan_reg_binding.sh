#!/bin/bash
echo STARTING PRC1 COOP SCAN
PARAM_FILE="params_static_overlap.yaml"
echo BASE PARAM FILE is $PARAM_FILE
for CONC_EFF in {10,15,20}; do
	for K_ON_SCALE in {500,525,550,575}; do
		K_ON=$(echo "scale=3; $K_ON_SCALE * 0.00001" | bc)
		FILE_NAME="static_overlapL"
		FILE_NAME+="_"
		FILE_NAME+=$K_ON_SCALE
		FILE_NAME+="_"
		FILE_NAME+=$CONC_EFF
		echo RUNNING NEW SIM: filename is $FILE_NAME
		TEMP_PARAMS="params_temp_static_"
		TEMP_PARAMS+=$K_ON_SCALE
		TEMP_PARAMS+="_"
		TEMP_PARAMS+=$CONC_EFF
		TEMP_PARAMS+=".yaml"
		cp $PARAM_FILE $TEMP_PARAMS
		yq w -i $TEMP_PARAMS xlinks.k_on $K_ON
		yq w -i $TEMP_PARAMS xlinks.c_eff_bind $CONC_EFF
		# Run sim for these parameter values
		./sim $TEMP_PARAMS $FILE_NAME &
	done
done
wait
rm params_temp_static_*
echo END SCAN
