#!/bin/bash
INITIAL_SEED=32443532;
N_SEEDS=1;
I_SEED=0;

while [ $I_SEED -lt $N_SEEDS ]; do
	SEED=$(echo "scale=0; $INITIAL_SEED + $I_SEED * 10" | bc)
	for MT_LENGTH in {250,500,750,1000,1250,1750}; do
		FILE_NAME="params_endtag_"
		FILE_NAME+=$I_SEED
		FILE_NAME+="_"
		FILE_NAME+=$MT_LENGTH
		FILE_NAME+=".yaml"
		cp params_base_endtag.yaml $FILE_NAME
		yq w -i $FILE_NAME microtubules.length[0] $MT_LENGTH
		yq w -i $FILE_NAME seed $SEED
		echo "created $FILE_NAME w/ MT = $MT_LENGTH sites"
	done
	let I_SEED+=1;
done
