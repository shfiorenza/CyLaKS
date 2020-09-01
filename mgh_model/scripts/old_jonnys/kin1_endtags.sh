#!/bin/bash
echo START KINESIN-1 ENDTAG
cd ~
cd Projects/overlap_analysis/mgh_model
#for i in {2,10,};
#do
#	XLINK_CONC=$(echo "scale=3; $i * 0.1" | bc)
	XLINK_CONC=1;
	for j in {8,12};
	do
		MOT_CONC=$(echo "scale=3; $j * 1.0" | bc)
		for k in {250,500,750,1000,1250,1750};
		do 
			FILE_NAME="Endtag"
			FILE_NAME+="_"
#			if [ $i -le 9 ]; then
#				FILE_NAME+="0"
#			fi
			FILE_NAME+="$XLINK_CONC"
			FILE_NAME+=".0"
			FILE_NAME+="x_"
			FILE_NAME+="$MOT_CONC"
			FILE_NAME+=".0"
			FILE_NAME+="m_"
			FILE_NAME+="$k"
			echo RUNNING NEW SIM: filename is $FILE_NAME
			# Update YAML files
			yq w -i params_endtags.yaml xlinks.concentration $XLINK_CONC
			yq w -i params_endtags.yaml motors.concentration $MOT_CONC
			yq w -i params_endtags.yaml microtubules.length $k
			# Run sim for these parameter values
			./sim params_endtags.yaml $FILE_NAME
			# Remove unnecessary files after each sim (save mem.)
			MT_COORD_FILE="$FILE_NAME"
			MT_COORD_FILE+="_mt_coord.file"
			rm $MT_COORD_FILE
			TETH_COORD_FILE="$FILE_NAME"
			TETH_COORD_FILE+="_tether_coord.file"
			rm $TETH_COORD_FILE
			MOTOR_FORCE_FILE="$FILE_NAME"
			MOTOR_FORCE_FILE+="_motor_force.file"
			rm $MOTOR_FORCE_FILE
			XLINK_FORCE_FILE="$FILE_NAME"
			XLINK_FORCE_FILE+="_xlink_force.file"
			rm $XLINK_FORCE_FILE
			TOT_FORCE_FILE="$FILE_NAME"
			TOT_FORCE_FILE+="_total_force.file"
			rm $TOT_FORCE_FILE
			MOTOR_EXT_FILE="$FILE_NAME"
			MOTOR_EXT_FILE+="_motor_extension.file"
			rm $MOTOR_EXT_FILE
			XLINK_EXT_FILE="$FILE_NAME"
			XLINK_EXT_FILE+="_xlink_extension.file"
			rm $XLINK_EXT_FILE
		done
	done
#done
echo END SCAN
