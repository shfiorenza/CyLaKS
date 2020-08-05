#!/bin/bash 
BASE_NAME="endtag_coop_full"
BASE_PARAMS="params_endtag.yaml"
echo "Starting ${BASE_NAME} scan"
echo "Base parameter file is ${BASE_PARAMS}"

BASE_SEED=198261346419
#for E_INT in 0 2 4 6 8 10
#do
	for N_SITES in 250 500 750 1000 1250 1750  
	do
		for SEED_NO in 0 1 2 3
		do
			SEED=$(( $BASE_SEED + $SEED_NO ))
	#		SIM_NAME="${BASE_NAME}_${E_INT}_${N_SITES}_${SEED_NO}"
			SIM_NAME="${BASE_NAME}_${N_SITES}_${SEED_NO}"
       			PARAM_FILE="params_temp_${SIM_NAME}.yaml"
        		echo "Launching sim ${SIM_NAME} with parameter file ${PARAM_FILE}"
      			cp ${BASE_PARAMS} ${PARAM_FILE}
    		    	yq w -i ${PARAM_FILE} seed ${SEED}
			yq w -i ${PARAM_FILE} microtubules.length[0] ${N_SITES}
		#	yq w -i ${PARAM_FILE} motors.interaction_energy ${E_INT}
	       	        # Run simulation; '&' allows for all to run concurrently 
     	   		./sim ${PARAM_FILE} ${SIM_NAME} & 
		done
	done
	wait
	rm params_temp_${BASE_NAME}_*
#done
