#!/bin/bash
INITIAL_SEED=32443532;
N_SEEDS=1;
I_SEED=0;

while [ $I_SEED -lt $N_SEEDS ]; do
	SEED=$(echo "scale=0; $INITIAL_SEED + $I_SEED * 10" | bc)
	for INITIAL_OVERLAP in {50,100,200,400,600,800,1000}; do
		FILE_NAME="params_slide_"
		FILE_NAME+=$I_SEED
		FILE_NAME+="_"
		FILE_NAME+=$INITIAL_OVERLAP
		FILE_NAME+=".yaml"
		cp params_base.yaml $FILE_NAME
		yq w -i $FILE_NAME microtubules.length[1] $INITIAL_OVERLAP
		yq w -i $FILE_NAME seed $SEED
		echo "created $FILE_NAME w/ top MT = $INITIAL_OVERLAP sites"
	done
	let I_SEED+=1;
done
