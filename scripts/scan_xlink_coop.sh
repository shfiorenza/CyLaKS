#!/bin/bash
BASE_NAME="prc1_coop"
BASE_PARAMS="params_endtag.yaml"
echo "Starting ${BASE_NAME} scan"
echo "Base parameter file is ${BASE_PARAMS}"

MIN_CONC=2
MAX_CONC=40
BASE_SEED=198261346419

 #for XLINK_CONC in 0 1 2 9 19 28 38
for (( XLINK_CONC=MIN_CONC; XLINK_CONC <= MAX_CONC; XLINK_CONC+=2 ))
do
	for I_SEED in 0 1 2 3 4 5 6 7 8 9
	do
		# SIM_NAME="${BASE_NAME}_${XLINK_CONC}"
		SIM_NAME="${BASE_NAME}_${XLINK_CONC}_${I_SEED}"
		PARAM_FILE="params_temp_${SIM_NAME}.yaml"
		echo "Launching sim ${SIM_NAME} with parameter file ${PARAM_FILE}"
		cp ${BASE_PARAMS} ${PARAM_FILE}
		yq w -i ${PARAM_FILE} seed $(( ${BASE_SEED} + ${I_SEED} ))
		yq w -i ${PARAM_FILE} xlinks.c_bulk $XLINK_CONC
		# Run sim for these parameter values
		./sim ${PARAM_FILE} ${SIM_NAME} &
	done
	wait
	rm params_temp_${BASE_NAME}_*
done
