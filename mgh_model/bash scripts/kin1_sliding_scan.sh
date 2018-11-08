#!/bin/bash
echo START KINESIN-1 SLIDING SCAN
cd ~
cd Projects/overlap_analysis/mgh_model
for k in `seq 16 20`;
do 
	FILE_NAME="SlideScan"
	FILE_NAME+="_"
	XLINK_CONC=$(echo "scale=2; 0.2 * $k" | bc)
	if [ $k -le 9 ]; then
		FILE_NAME+="0"
	fi
	FILE_NAME+="$XLINK_CONC"
	echo $FILE_NAME
	echo RUNNING NEW SIM: filename is $FILE_NAME
	# Update YAML files
	yq w -i params.yaml xlinks.concentration $XLINK_CONC
	# Run sim for these parameter values
	./sim params.yaml $FILE_NAME
	# Remove unnecessary files after each sim (save mem.)
	MOTOR_FORCE_FILE="$FILE_NAME"
	MOTOR_FORCE_FILE+="_motor_force.file"
	rm $MOTOR_FORCE_FILE
	MOTOR_EXT_FILE="$FILE_NAME"
	MOTOR_EXT_FILE+="_motor_extension.file"
	rm $MOTOR_EXT_FILE
	TOT_FORCE_FILE="$FILE_NAME"
	TOT_FORCE_FILE+="_total_force.file"
	rm $TOT_FORCE_FILE
done
echo END SCAN
