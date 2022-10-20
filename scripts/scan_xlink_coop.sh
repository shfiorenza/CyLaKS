#!/bin/bash
# BASE_NAME="xlink_coop"
BASE_NAME="shep_baseline"
BASE_PARAMS="params/k401.yaml"
echo "Starting ${BASE_NAME} scan"
echo "Base parameter file is ${BASE_PARAMS}"

MIN_CONC=2
MAX_CONC=40
BASE_SEED=198261346419

BASE_KON=0.000238
BASE_KOFF=0.143

#for (( XLINK_CONC=MIN_CONC; XLINK_CONC <= MAX_CONC; XLINK_CONC+=2 ))
#for XLINK_CONC in 15 30 45 60 75 90 105
for XLINK_AFF in 2 3 4 5 6 7 8 9 10
do
	for I_SEED in 0 #1 2 3 4 5 6 7 8 9
	do
		# SIM_NAME="${BASE_NAME}_${XLINK_CONC}nM"
		SIM_NAME="${BASE_NAME}_1nM_${XLINK_AFF}x"
		# SIM_NAME="${BASE_NAME}_${XLINK_CONC}_${I_SEED}"
		PARAM_FILE="params_temp_${SIM_NAME}.yaml"
		cp ${BASE_PARAMS} ${PARAM_FILE}
		KON=$(echo "scale=2; ${BASE_KON} * ${XLINK_AFF}" | bc)
		KOFF=$(echo "scale=2; ${BASE_KOFF} / ${XLINK_AFF}" | bc)
		echo "Launching sim ${SIM_NAME} with parameter file ${PARAM_FILE}"
    	yq eval -i ".seed = $(( ${BASE_SEED} + ${I_SEED} ))" ${PARAM_FILE}
	    # yq eval -i ".xlinks.c_bulk = ${XLINK_CONC}" ${PARAM_FILE}
	    yq eval -i ".xlinks.k_on = ${KON}" ${PARAM_FILE}
	    yq eval -i ".xlinks.k_off_i = ${KOFF}" ${PARAM_FILE}
		# Run sim for these parameter values
		./cylaks.exe ${PARAM_FILE} ${SIM_NAME} &
	done
	wait
done
rm params_temp_${BASE_NAME}_*
