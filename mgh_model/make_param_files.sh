#!/bin/bash
MT_LENGTH=1000;
STEP_SIZE=200;
N_FILES=$(echo "scale=0; (2 * $MT_LENGTH / $STEP_SIZE) + 1" | bc)
echo Making $N_FILES param files w/ $MT_LENGTH sites per MT
cd ~
cd Projects/overlap_analysis/mgh_model
yq w -i params_base.yaml microtubules.length $MT_LENGTH

COORD_SHIFT=0
while [	$COORD_SHIFT -le $MT_LENGTH ]; do
	FILE_NAME_0="params_shift_0_$COORD_SHIFT.yaml"
	FILE_NAME_1="params_shift_1_$COORD_SHIFT.yaml"
	cp params_base.yaml $FILE_NAME_0
	cp params_base.yaml $FILE_NAME_1
	yq w -i $FILE_NAME_0 microtubules.start_coord[0] $COORD_SHIFT
	yq w -i $FILE_NAME_1 microtubules.start_coord[1] $COORD_SHIFT
	echo "created $FILE_NAME_0 & $FILE_NAME_1 w/ shift=$COORD_SHIFT"
	let COORD_SHIFT+=$(STEP_SIZE);
done
