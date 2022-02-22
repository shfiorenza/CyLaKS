#!/bin/bash
BASE_NAME="k401_forceVel"
BASE_PARAMS="params/k401.yaml"
echo "Starting ${BASE_NAME} scan"
echo "Base parameter file is ${BASE_PARAMS}"

BASE_SEED=198261346419
for FORCE in 0 1 2 3 4 5 6
do
	for I_SEED in 0 1 2 3 4 5
	do
		SIM_NAME="${BASE_NAME}_${FORCE}_${I_SEED}"
		# SIM_NAME="${BASE_NAME}_${FORCE}"
		PARAM_FILE="params_temp_${SIM_NAME}.yaml"
		echo "Launching sim ${SIM_NAME} w/ parameter file ${PARAM_FILE}"
		cp ${BASE_PARAMS} ${PARAM_FILE}
	    yq eval -i ".seed = $(( ${BASE_SEED} + ${I_SEED} ))" ${PARAM_FILE}
	    yq eval -i ".motors.applied_force = ${FORCE}" ${PARAM_FILE}
		# Run sim for these parameter values
		./cylaks.exe ${PARAM_FILE} ${SIM_NAME} &
	done
done
wait
rm params_temp_${BASE_NAME}_*
