#!/bin/bash 
BASE_NAME="endtag_fluorescence"
BASE_PARAMS="params_endtag.yaml"
echo "Starting ${BASE_NAME} scan"
echo "Base parameter file is ${BASE_PARAMS}"

MOT_CONC=(0.02 1 2 4 6)
N_RUNS=(50 500 1000 2000 3000)

N_SIM_RUNS=0
N_CORES=12

#BASE_SEED=198261346419
#for E_INT in 0 2 4 6 8 10
for I_CONC in 0 1 2 3 4
do
	for N_SITES in 1250 # 250 500 750 1000 1250 1750  
	do
		#for SEED_NO in 0 1 2 3
		#do
	#		SEED=$(( $BASE_SEED + $SEED_NO ))
	#		SIM_NAME="${BASE_NAME}_${E_INT}_${N_SITES}_${SEED_NO}"
	#		SIM_NAME="${BASE_NAME}_${N_SITES}_${SEED_NO}"
			SIM_NAME="${BASE_NAME}_${MOT_CONC[I_CONC]}"
       		PARAM_FILE="params_temp_${SIM_NAME}.yaml"
        	echo "Launching sim ${SIM_NAME} with parameter file ${PARAM_FILE}"
      		cp ${BASE_PARAMS} ${PARAM_FILE}
    	#	yq w -i ${PARAM_FILE} seed ${SEED}
			yq w -i ${PARAM_FILE} microtubules.length[0] ${N_SITES}
			yq w -i ${PARAM_FILE} motors.c_bulk ${MOT_CONC[I_CONC]}
			yq w -i ${PARAM_FILE} motors.n_runs_desired ${N_RUNS[I_CONC]}
		#	yq w -i ${PARAM_FILE} motors.interaction_energy ${E_INT}
	       	# Run simulation; '&' allows for all to run concurrently 
     	   	./sim ${PARAM_FILE} ${SIM_NAME} & 
			N_SIM_RUNS=$(($N_SIM_RUNS + 1))
			echo "N_SIM_RUNS = ${N_SIM_RUNS}"
			if [ ${N_SIM_RUNS} -ge ${N_CORES} ]
			then
				echo "WAITING ..."
				wait
				N_SIM_RUNS=0;
			fi
		#done
	done
done
wait
rm params_temp_${BASE_NAME}_*
