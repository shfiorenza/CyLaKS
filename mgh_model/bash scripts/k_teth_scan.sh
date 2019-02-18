#!/bin/bash
echo START K_TETHER SCAN
cd ~
cd Projects/overlap_analysis/mgh_model
for k in {1,10,100,1000,10000};
#for k in {3,6,9,12,15,18};
do 
	FILE_NAME="TethScan"
	FILE_NAME+="_"
	K_TETH=$(echo "scale=2; 0.005 * $k" | bc)
	if [ $k -le 199 ]; then
		FILE_NAME+="0"
	fi
	FILE_NAME+="$K_TETH"
	echo RUNNING NEW SIM: filename is $FILE_NAME
	# Update YAML files
	yq w -i params_tethscan.yaml motors.k_tether $K_TETH
	# Run sim for these parameter values
	./sim params_tethscan.yaml $FILE_NAME
done
echo END SCAN
