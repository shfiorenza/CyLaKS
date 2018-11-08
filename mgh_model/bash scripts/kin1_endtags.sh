#!/bin/bash
echo START KINESIN-1 ENDTAG
cd ~
cd Projects/overlap_analysis/mgh_model
for i in {250,500,750,1000,1250,1750};
do 
	FILE_NAME="Endtag"
	FILE_NAME+="_"
	FILE_NAME+="$i"
	echo RUNNING NEW SIM: filename is $FILE_NAME
	# Update YAML files
	yq w -i params_endtags.yaml microtubules.length $i
	# Run sim for these parameter values
	./sim params_endtags.yaml $FILE_NAME
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
