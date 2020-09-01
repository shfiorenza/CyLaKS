#!/bin/bash
echo START MOTOR PROCESSIVITY SCAN
cd ~
cd Projects/overlap_analysis/mgh_model
for k in {0,1,2,3,4,5,6};
do 
	FILE_NAME="MotorAppForce"
	FILE_NAME+="_"
	APP_FORCE=$(echo "scale=2; 1.0 * $k" | bc)
	if [ $k -le 0 ]; then
		FILE_NAME+="0."
	fi
	FILE_NAME+="$APP_FORCE"
	FILE_NAME+="pN_b"
	echo RUNNING NEW SIM: filename is $FILE_NAME
	# Update YAML files
	yq w -i params.yaml motors.applied_force $APP_FORCE
	# Run sim for these parameter values
	./sim params.yaml $FILE_NAME
done
echo END SCAN
