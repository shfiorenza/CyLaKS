#!/bin/bash
BASE_NAME="kif4a_mobility"
BASE_PARAMS="params_processivity.yaml"
echo "Starting ${BASE_NAME} scan"
echo "Base parameter file is ${BASE_PARAMS}"

BASE_SEED=198261346419
for I_SEED in 0 1 2 3 4 5
do
	SIM_NAME="${BASE_NAME}_${I_SEED}"
	PARAM_FILE="params_temp_${SIM_NAME}.yaml"
	echo "Launching sim ${SIM_NAME} w/ parameter file ${PARAM_FILE}"
	cp ${BASE_PARAMS} ${PARAM_FILE}
	yq w -i ${PARAM_FILE} seed $(( ${BASE_SEED} + ${I_SEED} ))
	# Run sim for these parameter values
	./sim ${PARAM_FILE} ${SIM_NAME} &
done
wait
rm params_temp_${BASE_NAME}_*
echo END SCAN
