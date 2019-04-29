#!/bin/bash
INITIAL_SEED=32443532;
I_SEED=0;
N_SEEDS=24;
MT_LENGTH=1000;
STEP_SIZE=200;
N_FILES=$(echo "scale=0; $N_SEEDS * $MT_LENGTH / $STEP_SIZE" | bc)
echo Making $N_FILES param files w/ $MT_LENGTH sites per MT
yq w -i params_base.yaml microtubules.length $MT_LENGTH

while [ $I_SEED -lt $N_SEEDS ]; do
#	COORD_SHIFT=950
	SEED=$(echo "scale=0; $INITIAL_SEED + $I_SEED * 10" | bc)
#	while [	$COORD_SHIFT -lt $MT_LENGTH ]; do
	for COORD_SHIFT in {950,900,800,600,400,200,0}; do
		FILE_NAME="params_shift_stall_"
		FILE_NAME+=$I_SEED
		FILE_NAME+="_"
		FILE_NAME+=$COORD_SHIFT
		FILE_NAME+=".yaml"
		cp params_stall_base.yaml $FILE_NAME
		yq w -i $FILE_NAME microtubules.start_coord[0] $COORD_SHIFT
		yq w -i $FILE_NAME seed $SEED
		echo "created $FILE_NAME w/ shift=$COORD_SHIFT"
#		let COORD_SHIFT+=$STEP_SIZE;
	done
	let I_SEED+=1;
done
