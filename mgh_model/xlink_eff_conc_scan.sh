#!/bin/bash
echo START KINESIN-1 SLIDING SCAN
cd ~
cd Projects/overlap_analysis/mgh_model
# Run through desired MT lengths (2um, 4um, 6um, 8, 10um)
for k in `seq 1 5`;
do 
	FILE_NAME="XlinkDiffScan"
	FILE_NAME+="_"
	POW=1
	for i in seq `1 $k`;
	do
		POW=10*$POW
	done
	XLINK_CONC=$(echo "scale=2; 0.2 * $k" | bc)
	if [ $k -le 5 ]; then
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
	MOTOR_FILE="$FILE_NAME"
	MOTOR_FILE+="_motorID.file"
#	rm $MOTOR_FILE
	XLINK_FILE="$FILE_NAME"
	XLINK_FILE+="_xlinkID.file"
#	rm $XLINK_FILE
	MT_FILE="$FILE_NAME"
	MT_FILE+="_MTcoord.file"
#	rm $MT_FILE
done
echo END PARAMETER SCAN
