#!/bin/bash
BASE_NAME="hetero_tubulin"
BASE_PARAMS="params_processivity.yaml"
echo "Starting ${BASE_NAME} scan"
echo "Base parameter file is ${BASE_PARAMS}"

BASE_SEED=198261346419
for FRAC in 0 0.1 0.2 0.4 0.6 0.8 0.9 1
do
	for I_SEED in 0 1 2 3 4 5
	do
		SIM_NAME="${BASE_NAME}_${FRAC}_${I_SEED}"
		PARAM_FILE="params_temp_${SIM_NAME}.yaml"
		echo "Launching sim ${SIM_NAME} w/ parameter file ${PARAM_FILE}"
		cp ${BASE_PARAMS} ${PARAM_FILE}
		yq w -i ${PARAM_FILE} seed $(( ${BASE_SEED} + ${I_SEED} ))
		# Run sim for these parameter values
		./sim ${PARAM_FILE} ${SIM_NAME} ${FRAC} 2 &
	done
done
wait
rm params_temp_${BASE_NAME}_*
echo END SCAN
