#!/bin/bash
SCAN_NAME="lattice_coop"
echo Starting ${SCAN_NAME} scan
PARAM_FILE="params_processivity.yaml"
echo Base parameter file is ${PARAM_FILE}

N_SITES=(75000 50000 40000 25000 8500 3250)
MOT_CONC=(0.02 0.05 0.08 0.120 0.220 0.420)
CONC_SCALE=(20 50 80 120 220 420)

#for EMAX_BULK in 4.5
#do
#	for EMAX_SOLO in 0.0 1.5 2.0 2.5
#	do
	#	for CONC_SCALE in 20 50 80 120 220 420
			for I_CONC in 0 1 2 3 4 5
		do
		#	FILE_NAME="${SCAN_NAME}_${EMAX_BULK}_${EMAX_SOLO}_${CONC_SCALE}"
			FILE_NAME="${SCAN_NAME}_${CONC_SCALE[I_CONC]}"
			TEMP_PARAMS="params_temp_${FILE_NAME}.yaml"
			cp $PARAM_FILE $TEMP_PARAMS
			#N_STEPS=$(( 420 * 3600000 / $CONC_SCALE ))
#			MOT_CONC=$(echo "scale=4; $CONC_SCALE * 0.001" | bc)
			#yq w -i $TEMP_PARAMS n_steps $N_STEPS
			yq w -i $TEMP_PARAMS motors.c_bulk ${MOT_CONC[I_CONC]}
			yq w -i $TEMP_PARAMS microtubules.length[0] ${N_SITES[I_CONC]}
	#		yq w -i $TEMP_PARAMS motors.lattice_coop_Emax_solo $EMAX_SOLO
	#		yq w -i $TEMP_PARAMS motors.lattice_coop_Emax_bulk $EMAX_BULK
			# Run sim for these parameter values
			echo Running new sim: file name is ${FILE_NAME}
			./sim $TEMP_PARAMS $FILE_NAME &
		done
#	done
#done
wait
rm params_temp_${SCAN_NAME}*
echo END SCAN
