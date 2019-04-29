#!/bin/bash
echo STARTING PROCESSIVITY SCAN
for K_HYDROLYZE in {50,170}; do
	yq w -i params_processivity.yaml motors.k_hydrolyze $K_HYDROLYZE
for K_OFF_I in {30,6,3,2}; do
	C_BULK=$(echo "scale=5; 0.01 * $K_OFF_I / 30" | bc)
	PROCESSIVITY=$(echo "scale=1; 1.2 * 30 / $K_OFF_I" | bc)
	yq w -i params_processivity.yaml motors.c_bulk $C_BULK
	yq w -i params_processivity.yaml motors.k_off_i $K_OFF_I
#	for N_SITES in {250,500,750,1000,1250,1750}; do 
		FILE_NAME="Processivity"
		FILE_NAME+="_"
		FILE_NAME+="$PROCESSIVITY"
		FILE_NAME+="_"
		FILE_NAME+="$K_HYDROLYZE"
		echo RUNNING NEW SIM: filename is $FILE_NAME
		# Run sim for these parameter values
		./sim params_processivity.yaml $FILE_NAME
#		done
done
done
echo END SCAN
