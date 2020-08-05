#!/bin/bash
SCAN_NAME="kif4a_mobility_coop_full"
echo Starting ${SCAN_NAME} scan
PARAM_FILE="params_processivity.yaml"
echo Base parameter file is ${PARAM_FILE}

MOT_CONC=(0.02 0.05 0.08 0.120 0.220 0.420)
CONC_SCALE=(20 50 80 120 220 420)
#N_RUNS=(7 18 30 50 100 150)
N_RUNS=(50 125 200 300 500 800)
#N_RUNS=(500 1250 2000 3000 5000 8000)

BASE_SEED=198261346419

#for E_INT in 2 4 5 6 8 10
#do
	for I_CONC in 0 1 2 3 4 5
	do
		for SEED_NO in 0 1 2 3
		do
			SEED=$(( $BASE_SEED + $SEED_NO ))
		#	FILE_NAME="${SCAN_NAME}_${E_INT}_${CONC_SCALE[I_CONC]}_${SEED_NO}"
			FILE_NAME="${SCAN_NAME}_${CONC_SCALE[I_CONC]}_${SEED_NO}"
		#	FILE_NAME="${SCAN_NAME}_${CONC_SCALE[I_CONC]}"
			TEMP_PARAMS="params_temp_${FILE_NAME}.yaml"
			cp $PARAM_FILE $TEMP_PARAMS
			yq w -i ${TEMP_PARAMS} seed ${SEED}
			yq w -i ${TEMP_PARAMS} motors.c_bulk ${MOT_CONC[I_CONC]}
			yq w -i ${TEMP_PARAMS} motors.n_runs_desired ${N_RUNS[I_CONC]}
		#	yq w -i ${TEMP_PARAMS} motors.interaction_energy ${E_INT}
			# Run sim for these parameter values
			echo Running new sim: file name is ${FILE_NAME}
			./sim $TEMP_PARAMS $FILE_NAME &
		done
	done
	wait
	rm params_temp_${SCAN_NAME}*
#done
echo END SCAN
