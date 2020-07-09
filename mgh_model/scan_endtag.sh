#!/bin/bash 
BASE_NAME="endtag"
BASE_PARAMS="params_processivity.yaml"
echo "Starting ${BASE_NAME} scan"
echo "Base parameter file is ${BASE_PARAMS}"

BASE_SEED=198261346419
for CONC_SCALE in 0
do
	XLINK_CONC=$(echo "scale=2; $CONC_SCALE * 0.1" | bc)
	for N_SITES in 250 500 750 1000 1250 1750  
	do
		for SEED_NO in 0 1 2 3 4 5 6 7 8 9 
		do
			SEED=$(( $BASE_SEED + $SEED_NO ))
			SIM_NAME="${BASE_NAME}_${N_SITES}_${SEED_NO}"
       		PARAM_FILE="params_temp_${SIM_NAME}.yaml"
        	echo "Launching sim ${SIM_NAME} with parameter file ${PARAM_FILE}"
      		cp ${BASE_PARAMS} ${PARAM_FILE}
    	    yq w -i ${PARAM_FILE} seed ${SEED}
			yq w -i ${PARAM_FILE} xlinks.c_bulk ${XLINK_CONC}
			yq w -i ${PARAM_FILE} microtubules.length[0] ${N_SITES}
     	   # Run simulation; '&' allows for all to run concurrently 
     	   ./sim ${PARAM_FILE} ${SIM_NAME} & 
		done
		wait
		rm params_temp_${BASE_NAME}_*
	done
done
