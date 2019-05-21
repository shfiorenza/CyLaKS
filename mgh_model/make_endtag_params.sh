#!/bin/bash
#INITIAL_SEED=32443532;
#N_SEEDS=1;
#I_SEED=0;

#while [ $I_SEED -lt $N_SEEDS ]; do
#	SEED=$(echo "scale=0; $INITIAL_SEED + $I_SEED * 10" | bc)
#for K_HYDRO in {75,90,105,120,135}; do
K_HYDRO=60;
	for JAM_RATIO in {100,200,300,400,500,600,700,800}; do
		K_HYDRO_ST=$(echo "scale=3; $K_HYDRO / $JAM_RATIO" | bc)
		for MT_LENGTH in {250,500,750,1000,1250,1750}; do
			FILE_NAME="params_endtag_"
			FILE_NAME+=$K_HYDRO
			FILE_NAME+="_"
			FILE_NAME+=$JAM_RATIO
			FILE_NAME+="_"
			FILE_NAME+=$MT_LENGTH
			FILE_NAME+=".yaml"
			cp params_base_endtag.yaml $FILE_NAME
			yq w -i $FILE_NAME motors.k_hydrolyze $K_HYDRO
			yq w -i $FILE_NAME motors.k_hydrolyze_stalled $K_HYDRO_ST
			yq w -i $FILE_NAME microtubules.length[0] $MT_LENGTH
#			yq w -i $FILE_NAME seed $SEED
			echo "created param file '$FILE_NAME'"
		done
	done
#	let I_SEED+=1;
#done
