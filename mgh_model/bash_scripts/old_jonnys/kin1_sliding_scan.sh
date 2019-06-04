#!/bin/bash
echo START KINESIN-1 SLIDING SCAN
cd ~
cd Projects/overlap_analysis/mgh_model
for k in {1,10,100,1000};
do 
	FILE_NAME="SlideScan"
	FILE_NAME+="_"
	K_TETH=$(echo "scale=2; 0.01 * $k" | bc)
	if [ $k -le 99 ]; then
		FILE_NAME+="0"
	fi
	FILE_NAME+="$K_TETH"
	echo $FILE_NAME
	echo RUNNING NEW SIM: filename is $FILE_NAME
	# Update YAML files
	yq w -i params_slidescan.yaml motors.k_tether_free $K_TETH
	# Run sim for these parameter values
	./sim params_slidescan.yaml $FILE_NAME
done
echo END SCAN
