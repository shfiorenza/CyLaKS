#!/bin/bash
SCAN_NAME="endtag"
echo Starting ${SCAN_NAME} scan
PARAM_FILE="params_endtag.yaml"
echo Base parameter file is ${PARAM_FILE}

BASE_SEED=198261346419

for N_SITES in 250 500 750 1000 1250 1750
do
	for RANGE in 10 50 100 500 1000
	do
		for SEED_NO in 0 1 2 3
		do
			FILE_NAME="${SCAN_NAME}_${N_SITES}_${RANGE}_${SEED_NO}"
			TEMP_PARAMS="params_temp_${FILE_NAME}.yaml"
			cp $PARAM_FILE $TEMP_PARAMS
			yq w -i ${TEMP_PARAMS} seed $(( $BASE_SEED + $SEED_NO ))
			yq w -i ${TEMP_PARAMS} motors.gaussian_range ${RANGE}
			yq w -i ${TEMP_PARAMS} filaments.n_sites[0] ${N_SITES}
			# Run sim for these parameter values
			echo Running new sim: file name is ${FILE_NAME}
			./sim $TEMP_PARAMS $FILE_NAME &
		done
	done
done
wait
rm params_temp_${SCAN_NAME}*
echo END SCAN
