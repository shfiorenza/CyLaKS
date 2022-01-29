#!/bin/bash
BASE_NAME="separation"
BASE_PARAMS="params/overlap.yaml"
echo "Starting ${BASE_NAME} scan"
echo "Base parameter file is ${BASE_PARAMS}"

MIN_NUM=1
MAX_NUM=13
N_SEEDS=5

BASE_SEED=198261346419
for N_XLINKS in $(seq ${MIN_NUM} ${MAX_NUM})
do
	for I_SEED in $(seq 0 $((${N_SEEDS}-1)))
	do
		SIM_NAME="${BASE_NAME}_${N_XLINKS}_${I_SEED}"
		PARAM_FILE="temp_params_${SIM_NAME}.yaml"
		echo "Launching sim ${SIM_NAME} with parameter file ${PARAM_FILE}"
		cp ${BASE_PARAMS} ${PARAM_FILE}
	    yq eval -i ".seed = $(( ${BASE_SEED} + ${I_SEED} ))" ${PARAM_FILE}
		# Run sim for these parameter values
		./cylaks.exe ${PARAM_FILE} ${SIM_NAME} filament_separation ${N_XLINKS}
		wait
		rm temp_params_${BASE_NAME}_*
	done
done
