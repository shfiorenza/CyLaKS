#!/bin/bash
BASE_NAME="motor_heterodimer"
BASE_PARAMS="params/kif4a.yaml"
echo "Starting ${BASE_NAME} scan"
echo "Base parameter file is ${BASE_PARAMS}"

BASE_SEED=198261346419
for D in 0.01 0.05 0.1 0.5 1 5
do
	for I_SEED in 0 1 2 3
	do
		SIM_NAME="${BASE_NAME}_${D}_${I_SEED}"
		PARAM_FILE="params_temp_${SIM_NAME}.yaml"
		echo "Launching sim ${SIM_NAME} w/ parameter file ${PARAM_FILE}"
		cp ${BASE_PARAMS} ${PARAM_FILE}
	    yq eval -i ".seed = $(( ${BASE_SEED} + ${I_SEED} ))" ${PARAM_FILE}
   		yq eval -i ".xlinks.d_i = ${D}" ${PARAM_FILE}
		# Run sim for these parameter values
		./cylaks.exe ${PARAM_FILE} ${SIM_NAME} kinesin_mutant &
	done
done
wait
rm params_temp_${BASE_NAME}_*
