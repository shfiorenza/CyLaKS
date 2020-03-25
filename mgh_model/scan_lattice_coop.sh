#!/bin/bash
SCAN_NAME="lattice_coop"
echo Starting ${SCAN_NAME} scan
PARAM_FILE="params_processivity.yaml"
echo Base parameter file is ${PARAM_FILE}

for EMAX_BULK in 4.5
do
	for EMAX_SOLO in 0.0 1.5 2.0 2.5
	do
		for CONC_SCALE in 20 # 220 420
		do
			FILE_NAME="${SCAN_NAME}_${EMAX_BULK}_${EMAX_SOLO}_${CONC_SCALE}"
			echo Running new sim: file name is ${FILE_NAME}
			TEMP_PARAMS="params_temp_${FILE_NAME}.yaml"
			cp $PARAM_FILE $TEMP_PARAMS
			#N_STEPS=$(( 420 * 3600000 / $CONC_SCALE ))
			MOT_CONC=$(echo "scale=4; $CONC_SCALE * 0.001" | bc)
			#yq w -i $TEMP_PARAMS n_steps $N_STEPS
			yq w -i $TEMP_PARAMS motors.c_bulk $MOT_CONC
			yq w -i $TEMP_PARAMS motors.lattice_coop_Emax_solo $EMAX_SOLO
			yq w -i $TEMP_PARAMS motors.lattice_coop_Emax_bulk $EMAX_BULK
			# Run sim for these parameter values
			./sim $TEMP_PARAMS $FILE_NAME &
		done
	done
done
wait
rm params_temp_${SCAN_NAME}*
echo END SCAN
