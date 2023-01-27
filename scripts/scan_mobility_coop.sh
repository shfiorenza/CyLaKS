#!/bin/bash
SCAN_NAME="kif4a_mobility"
PARAM_FILE="params/kif4a.yaml"
echo Starting ${SCAN_NAME} scan
echo Base parameter file is ${PARAM_FILE}

BASE_SEED=198261346419

MOT_CONC=(0.02 0.05 0.08 0.120 0.220 0.420)
CONC_SCALE=(20 50 80 120 220 420)
#N_RUNS=(50 125 200 300 500 800)
N_RUNS=(20 50 80 120 220 420)

for RANGE in 5 10 50 100
do
	for I_CONC in 0 1 2 3 4 5
	do
		for SEED_NO in 0  1 2 3
		do
			FILE_NAME="${SCAN_NAME}_${RANGE}_${CONC_SCALE[I_CONC]}_${SEED_NO}"
			# FILE_NAME="${SCAN_NAME}_${CONC_SCALE[I_CONC]}"
			TEMP_PARAMS="params_temp_${FILE_NAME}.yaml"
			cp $PARAM_FILE $TEMP_PARAMS
    		yq eval -i ".seed = $(( $BASE_SEED + $SEED_NO ))}" ${TEMP_PARAMS}
    		yq eval -i ".motors.n_runs_to_exit = ${N_RUNS[I_CONC]}" ${TEMP_PARAMS}
			yq eval -i ".motors.c_bulk = ${MOT_CONC[I_CONC]}" ${TEMP_PARAMS}
			yq eval -i ".motors.gaussian_range = ${RANGE}" ${TEMP_PARAMS}
			# Run sim for these parameter values
			echo Running new sim: file name is ${FILE_NAME}
			./cylaks.exe $TEMP_PARAMS $FILE_NAME &
		done
	done
	wait
	rm params_temp_${SCAN_NAME}*
done
