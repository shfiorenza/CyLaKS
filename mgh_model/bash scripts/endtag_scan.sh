#!/bin/bash
echo STARTING ENDTAG SCAN
#MT_LENGTH=1000;
#yq w -i params_endtags.yaml microtubules.length $MT_LENGTH
#for i in {5,10,20}; do
#	XLINK_CONC=$(echo "scale=3; $i * 0.1" | bc)
#	for j in {20,100,200}; do
#		MOT_CONC=$(echo "scale=3; $j * 1.0" | bc)
mkdir output_new
#for K_ON in {1,5,10,15}; do
#	K_ON=$(echo "scale=4; $K_ON / 100" | bc);
#	yq w -i params_endtags.yaml motors.k_on $K_ON
#	for K_OFF_I in {30,6,3,2,1}; do
#		PROCESSIVITY=$(echo "scale=1; 1.2 * 30 / $K_OFF_I" | bc)
#		C_BULK=$(echo "scale=4; 1.5 * $K_OFF_I / 30" | bc)
#		yq w -i params_endtags.yaml motors.k_off_i $K_OFF_I
#		yq w -i params_endtags.yaml motors.c_bulk $C_BULK
		for N_SITES in {250,500,750,1000,1250,1750}; do 
			FILE_NAME="Endtag"
			FILE_NAME+="_"
#			FILE_NAME+="$PROCESSIVITY"
#			FILE_NAME+="_0"
#			FILE_NAME+="$K_ON"
#			FILE_NAME+="_"
#			if [ $i -le 9 ]; then
#				FILE_NAME+="0"
#			fi
#			if [ $i -eq 0 ]; then
#				FILE_NAME+="."
#			fi
#			FILE_NAME+="$XLINK_CONC"
#			FILE_NAME+=".0"
#			FILE_NAME+="x_"
#			FILE_NAME+="$MOT_CONC"
#			FILE_NAME+=".0"
#			FILE_NAME+="m_"
			FILE_NAME+="$N_SITES"
			echo RUNNING NEW SIM: filename is $FILE_NAME
			# Update YAML files
#			yq w -i params_endtags.yaml xlinks.concentration $XLINK_CONC
#			yq w -i params_endtags.yaml motors.c_bulk $MOT_CONC
			yq w -i params_endtags.yaml microtubules.length[0] $N_SITES
			# Run sim for these parameter values
			./sim params_endtags.yaml $FILE_NAME
			mv $FILE_NAME* output_new/
#		done
#	done
done
echo END SCAN
