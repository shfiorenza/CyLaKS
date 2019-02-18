#!/bin/bash
echo START XLINK DIFFUSION SCAN
cd ~
cd Projects/overlap_analysis/mgh_model
for k in {1,5,10,50,100,500};
do 
	echo $k
	FILE_NAME="XlinkConcScan"
	FILE_NAME+="_"
	CONC=$(echo "scale=3; 0.01 * $k" | bc)
	if [ $k -le 99 ]; then
		FILE_NAME+="0"
	fi
	FILE_NAME+="$CONC"
	echo RUNNING NEW SIM: filename is $FILE_NAME
	# Update YAML files
	yq w -i xlscan_params.yaml xlinks.concentration $CONC
	# Run sim for these parameter values
	./sim xlscan_params.yaml $FILE_NAME
done
echo END SCAN
