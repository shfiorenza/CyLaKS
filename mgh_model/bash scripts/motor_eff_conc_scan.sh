#!/bin/bash
echo START MOTOR PROCESSIVITY SCAN
cd ~
cd Projects/overlap_analysis/mgh_model
for k in `seq 10 15`;
do 
	echo $k
	FILE_NAME="MotorProScan"
	FILE_NAME+="_"
	POW=1
	for i in `seq 1 $(echo "scale=2; $k / 2" | bc)`;
	do
		POW=$(echo "scale=2; 10 * $POW" | bc)
	done
	EFF_CONC=$(echo "scale=2; 80 * $POW" | bc)
	if [ $(($k % 2)) -eq 1 ]; then
		EFF_CONC=$(echo "scale=2; 2 * $EFF_CONC" | bc)
	fi
	FILE_NAME+="$EFF_CONC"
	echo $FILE_NAME
	echo RUNNING NEW SIM: filename is $FILE_NAME
	# Update YAML files
	yq w -i motscan_params.yaml motors.conc_eff_bind $EFF_CONC
	# Run sim for these parameter values
	./sim motscan_params.yaml $FILE_NAME
	# Remove unnecessary files after each sim (save mem.)
	MT_FILE="$FILE_NAME"
	MT_FILE+="_mt_coord.file"
	rm $MT_FILE
	XLINK_FILE="$FILE_NAME"
	XLINK_FILE+="_xlinkID.file"
	rm $XLINK_FILE
	XLINK_FORCE_FILE="$FILE_NAME"
	XLINK_FORCE_FILE+="_xlink_force.file"
	rm $XLINK_FORCE_FILE
	XLINK_EXT_FILE="$FILE_NAME"
	XLINK_EXT_FILE+="_xlink_extension.file"
	rm $MOTOR_EXT_FILE
	TETH_FILE="$FILE_NAME"
	TETH_FILE+="_tether_coord.file"
	rm $TETH_FILE
	TOT_FORCE_FILE="$FILE_NAME"
	TOT_FORCE_FILE+="_total_force.file"
	rm $TOT_FORCE_FILE

done
echo END SCAN
