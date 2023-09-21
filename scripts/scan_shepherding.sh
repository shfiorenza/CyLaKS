#!/bin/bash
BASE_NAME="shep"
BASE_PARAMS="params/k401.yaml"
# BASE_PARAMS="params/slide.yaml"
echo "Starting ${BASE_NAME} scan"
echo "Base parameter file is ${BASE_PARAMS}"

BASE_SEED=198261346419

for XLINK_CONC in 1 10
do
	for MOT_CONC in 0
	do
	for E_INT in 0.0 0.75 1.5  
	do
		for N_PFS in 1 2 4 8 
		do
		# for DIFF_CONST in 0.131 # 0.0131 0.00131 0.000131 0.0000131
		# do
	#	for OFF_RATE in 0.0143 # 0.00143
	# 	do
	#	for N_SITES in 1000 1250 1750			
	#	do
	#  for I_SEED in 0 1 2 # 3 # 0 1 # 2 # 3
			# SIM_NAME="${BASE_NAME}_${XLINK_CONC}nM_${MOT_CONC}nM_${DIFF_CONST}_${N_PFS}"
			# SIM_NAME="${BASE_NAME}_${DIFF_CONST}_${N_PFS}"
			# SIM_NAME="${BASE_NAME}_${XLINK_CONC}nM_${MOT_CONC}nM_${DIFF_CONST}_${OFF_RATE}_${N_SITES}_short"
			SIM_NAME="${BASE_NAME}_${XLINK_CONC}nM_${MOT_CONC}nM_${N_PFS}_${E_INT}kT"
			PARAM_FILE="params_temp_${SIM_NAME}.yaml"
			cp ${BASE_PARAMS} ${PARAM_FILE}
			yq eval -i ".xlinks.c_bulk = ${XLINK_CONC}" ${PARAM_FILE}
			yq eval -i ".motors.c_bulk = ${MOT_CONC}" ${PARAM_FILE}
			# yq eval -i ".xlinks.d_i = ${DIFF_CONST}" ${PARAM_FILE}
			# yq eval -i ".xlinks.d_side = ${DIFF_CONST}" ${PARAM_FILE}
			yq eval -i ".filaments.n_subfilaments = ${N_PFS}" ${PARAM_FILE}
			yq eval -i ".xlinks.neighb_neighb_energy = ${E_INT}" ${PARAM_FILE}
			# yq eval -i ".xlinks.k_off_i = ${OFF_RATE}" ${PARAM_FILE}
			# yq eval -i ".filaments.n_sites[0] = ${N_SITES}" ${PARAM_FILE}
			# yq eval -i ".filaments.n_sites[1] = ${N_SITES}" ${PARAM_FILE}
		    # yq eval -i ".seed = $(( ${BASE_SEED} + ${I_SEED} ))" ${PARAM_FILE}
			# Run sim for these parameter values
			echo "Launching sim ${SIM_NAME} with parameter file ${PARAM_FILE}"
			# ./cylaks.exe ${PARAM_FILE} ${SIM_NAME} &
			singularity exec --bind $PWD cylaks_latest.sif cylaks.exe ${PARAM_FILE} ${SIM_NAME} &

	    # done
	#   done
		done
		done
	done
done
wait
# rm params_temp_${BASE_NAME}_*