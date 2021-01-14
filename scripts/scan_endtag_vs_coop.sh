#!/bin/bash
SCAN_NAME="endtag"
echo Starting ${SCAN_NAME} scan
PARAM_FILE="params_endtag.yaml"
echo Base parameter file is ${PARAM_FILE}

#MOT_CONC=(0.02 0.05 0.08 0.120 0.220 0.420)
MOT_CONC=(1 2 4 6)
#CONC_SCALE=(20 50 80 120 220 420)
CONC_SCALE=(1 2 4 6)
BASE_SEED=198261346419

for I_CONC in 0 1 2 3 # 4 5
do
	for RANGE in 10 100 1000
	do
		for SEED_NO in 0 1 2 # 3
		do
			SEED=$(( $BASE_SEED + $SEED_NO ))
			FILE_NAME="${SCAN_NAME}_${CONC_SCALE[I_CONC]}_${RANGE}_${SEED_NO}"
			TEMP_PARAMS="params_temp_${FILE_NAME}.yaml"
			cp $PARAM_FILE $TEMP_PARAMS
			yq w -i ${TEMP_PARAMS} seed ${SEED}
			yq w -i ${TEMP_PARAMS} motors.gaussian_range ${RANGE}
			yq w -i ${TEMP_PARAMS} motors.c_bulk ${MOT_CONC[I_CONC]}
			# Run sim for these parameter values
			echo Running new sim: file name is ${FILE_NAME}
			./sim $TEMP_PARAMS $FILE_NAME &
		done
	done
	wait
	rm params_temp_${SCAN_NAME}*
done
echo END SCAN
