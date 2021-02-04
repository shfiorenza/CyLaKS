#!/bin/bash
BASE_NAME="mt_forceVel"
BASE_PARAMS="params_slide.yaml"
echo "Starting ${BASE_NAME} scan"
echo "Base parameter file is ${BASE_PARAMS}"

BASE_SEED=198261346419

for APPLIED_FORCE in 50.0 5.0 0.5 0.05
do
    for I_SEED in 0 # 1 2 3 4 5
    do
        SEED=$(( ${BASE_SEED} + ${I_SEED} ))
        SIM_NAME="${BASE_NAME}_${APPLIED_FORCE}_${I_SEED}"
        PARAM_FILE="temp_params_${SIM_NAME}.yaml"
        echo "Launching sim ${SIM_NAME} with parameter file ${PARAM_FILE}"
        cp ${BASE_PARAMS} ${PARAM_FILE}
        yq eval -i ".filaments.f_applied[0] = ${APPLIED_FORCE}" ${PARAM_FILE}
        yq eval -i ".filaments.f_applied[1] = ${APPLIED_FORCE}" ${PARAM_FILE}
        yq eval -i ".filaments.f_applied[2] = ${APPLIED_FORCE}" ${PARAM_FILE}
        yq eval -i ".seed = ${SEED}" ${PARAM_FILE}
        # Run simulation; '&' allows for all to run concurrently 
        ./sim ${PARAM_FILE} ${SIM_NAME} & 
    done
done
wait 
 rm temp_params_${BASE_NAME}_*