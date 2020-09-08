#!/bin/bash 
echo "*** Starting endtag scan script ***"
PARAM_FILE="params_endtag.yaml"
echo "	Base parameter file is $PARAM_FILE"

BASE_SEED=198261346419
for CONC_SCALE in 0
do
	XLINK_CONC=$(echo "scale=2; $CONC_SCALE * 0.1" | bc)
	for N_SITES in 250 500 750 1000 1250 1750  
	do
		for SEED_NO in 0 1 2 3 4 5 6 7 8 9 
		do
			SEED=$(( $BASE_SEED + $SEED_NO ))
			FILE_NAME="endtagE"
			FILE_NAME+="_"
			FILE_NAME+=$CONC_SCALE
			FILE_NAME+="_"
			FILE_NAME+=$N_SITES
			FILE_NAME+="_"
			FILE_NAME+=$SEED_NO
			echo "*** Running new sim; file name is $FILE_NAME ***"
			TEMP_PARAMS="params_temp_"
			TEMP_PARAMS+=$FILE_NAME
			TEMP_PARAMS+=".yaml"
			cp $PARAM_FILE $TEMP_PARAMS
			yq w -i $TEMP_PARAMS seed $SEED
			yq w -i $TEMP_PARAMS xlinks.c_bulk $XLINK_CONC
			yq w -i $TEMP_PARAMS microtubules.length[0] $N_SITES
			# Run sim for these parameter values
			./sim $TEMP_PARAMS $FILE_NAME &
		done
		wait
		rm params_temp_*
	done
done
echo "*** Endtag scan script has finished ***"
