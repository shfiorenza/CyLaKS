#!/bin/bash 
BASE_NAME="endtag_fluorescence"
BASE_PARAMS="params/endtag.yaml"
echo "Starting ${BASE_NAME} scan"
echo "Base parameter file is ${BASE_PARAMS}"

MOT_CONC=(0.02 1 2 4 6)
N_RUNS=(50 500 1000 2000 3000)

N_SIM_RUNS=0
N_CORES=12

BASE_SEED=198261346419

#for E_INT in 0 2 4 6 8 10
for I_CONC in 0 1 2 3 4
do
	for N_SITES in 1250 # 250 500 750 1000 1250 1750  
	do
		#for SEED_NO in 0 1 2 3
		#do
	#		SIM_NAME="${BASE_NAME}_${E_INT}_${N_SITES}_${SEED_NO}"
	#		SIM_NAME="${BASE_NAME}_${N_SITES}_${SEED_NO}"
			SIM_NAME="${BASE_NAME}_${MOT_CONC[I_CONC]}"
       		PARAM_FILE="params_temp_${SIM_NAME}.yaml"
        	echo "Launching sim ${SIM_NAME} with parameter file ${PARAM_FILE}"
      		cp ${BASE_PARAMS} ${PARAM_FILE}
		    # yq eval -i ".seed = $(( ${BASE_SEED} + ${I_SEED} ))" ${PARAM_FILE}
			yq eval -i ".filaments.n_sites[0] = ${N_SITES}" ${PARAM_FILE}
			yq eval -i ".motors.c_bulk = ${MOT_CONC[I_CONC]}" ${PARAM_FILE}
    		yq eval -i ".motors.n_runs_to_exit = ${N_RUNS[I_CONC]}" ${PARAM_FILE}
			yq eval -i ".motors.interaction_energy = ${E_INT}" ${PARAM_FILE} 
	       	# Run simulation; '&' allows for all to run concurrently 
     	   	./cylaks.exe ${PARAM_FILE} ${SIM_NAME} & 
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
