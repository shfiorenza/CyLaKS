#!/bin/bash
echo STARTING PRC1 COOP SLIDING SCAN
PARAM_FILE="params_overlap.yaml"
echo BASE PARAM FILE is $PARAM_FILE
for C_EFF_TETH in {1500,5000}; do
for E_SCALE in {0,75,150,225,300}; do #,250,275,300,325}; do
	E_INT=$(echo "scale=3; $E_SCALE * 0.01" | bc)
		FILE_NAME="coop_slide_"
		FILE_NAME+=$C_EFF_TETH
		FILE_NAME+="_"
		FILE_NAME+=$E_SCALE
		echo RUNNING NEW SIM: filename is $FILE_NAME
		TEMP_PARAMS="params_temp_"
		TEMP_PARAMS+=$C_EFF_TETH
		TEMP_PARAMS+="_"
		TEMP_PARAMS+=$E_SCALE
		TEMP_PARAMS+=".yaml"
		cp $PARAM_FILE $TEMP_PARAMS
		yq w -i $TEMP_PARAMS xlinks.interaction_energy $E_INT
		yq w -i $TEMP_PARAMS motors.c_eff_tether $C_EFF_TETH
		# Run sim for these parameter values
		./sim $TEMP_PARAMS $FILE_NAME &
done
done
wait
rm params_temp_*
echo END SCAN