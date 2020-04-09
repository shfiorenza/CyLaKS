#!/bin/bash
BASE_NAME="xlink_diffusion_double"
BASE_PARAMS="params_slide.yaml"
echo "Starting ${BASE_NAME} scan"
echo "Base parameter file is ${BASE_PARAMS}"

#BASE_SEED=198261346419

#for I_SEED in 0 1 2 3 4 5 6 7
for I_MT in 0 1
do
	for OFFSET in 0 0.25 0.5 0.75 1
	do
		# SEED=$(( ${BASE_SEED} + ${I_SEED} ))
		# SIM_NAME="${BASE_NAME}_${I_SEED}"
		SIM_NAME="${BASE_NAME}_${I_MT}_${OFFSET}"
		PARAM_FILE="temp_params_${SIM_NAME}.yaml"
		cp ${BASE_PARAMS} ${PARAM_FILE}
		# yq w -i ${PARAM_FILE} seed ${SEED}
		yq w -i ${PARAM_FILE} microtubules.start_coord[${I_MT}] ${OFFSET}
		# Run simulation; '&' allows for all to run concurrently 
		echo "Launching sim ${SIM_NAME} with parameter file ${PARAM_FILE}"
		./sim ${PARAM_FILE} ${SIM_NAME} & 
	done
done
wait 
rm temp_params_${BASE_NAME}_*
echo "Finished ${BASE_NAME} scan"