#!/bin/bash
echo START K_UNTETHER SCAN
cd ~
cd Projects/overlap_analysis/mgh_model
for k in {1,10,100,1000,10000};
do 
	FILE_NAME="UntethScan"
	FILE_NAME+="_"
	K_TETH=$(echo "scale=2; 0.005 * $k" | bc)
	if [ $k -le 199 ]; then
		FILE_NAME+="0"
	fi
	FILE_NAME+="$K_TETH"
	echo $FILE_NAME
	echo RUNNING NEW SIM: filename is $FILE_NAME
	# Update YAML files
	yq w -i params_tethscan2.yaml motors.k_untether $K_TETH
	# Run sim for these parameter values
	./sim params_tethscan2.yaml $FILE_NAME
done
echo END SCAN
