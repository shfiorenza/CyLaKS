#!/bin/bash
echo START KINESIN-1 ENDTAG
for i in {0,1,4}; do
	XLINK_CONC=$(echo "scale=3; $i * 0.1" | bc)
#	for j in {8,12}; do
#		MOT_CONC=$(echo "scale=3; $j * 1.0" | bc)
		for k in {250,500,750,1000,1250,1750}; do 
			FILE_NAME="Endtag"
			FILE_NAME+="_"
			if [ $i -le 9 ]; then
				FILE_NAME+="0"
			fi
			FILE_NAME+="$XLINK_CONC"
			FILE_NAME+=".0"
			FILE_NAME+="x_"
#			FILE_NAME+="$MOT_CONC"
#			FILE_NAME+=".0"
#			FILE_NAME+="m_"
			FILE_NAME+="$k"
			echo RUNNING NEW SIM: filename is $FILE_NAME
			# Update YAML files
			yq w -i params_endtags.yaml xlinks.concentration $XLINK_CONC
#			yq w -i params_endtags.yaml motors.concentration $MOT_CONC
			yq w -i params_endtags.yaml microtubules.length $k
			# Run sim for these parameter values
			./sim params_endtags.yaml $FILE_NAME
		done
#	done
done
echo END SCAN
