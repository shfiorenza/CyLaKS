#!/bin/bash
echo START PARAMETER SCAN
cd ~
cd Projects/overlap_analysis/mgh_model
# Run through different k_off ratios: 5, 25, 50, 250, 500, 2500
for h in `seq 6 6`;
do
	K_OFF_RATIO=5
	if [ $(($h % 2)) -eq 1 ]; then
		for pow in `seq 2 2 $h`;
		do
			K_OFF_RATIO=$(($K_OFF_RATIO * 10))
		done
	else
		K_OFF_RATIO=$(($K_OFF_RATIO * 5));
		for pow in `seq 3 2 $h`;
		do
			K_OFF_RATIO=$(($K_OFF_RATIO * 10))
		done
	fi
	# Run through different failstep rates: 0.005, 0.05, 0.5, 5, 25, 50
	for i in `seq 1 6`;
	do
		FAILSTEP_RATE=$(echo "scale=3 ; 0.005" | bc)
		if [ $i -lt 5 ]; then
			for pow in `seq 1 $(($i - 1))`;
			do
				FAILSTEP_RATE=$(echo "scale=3 ; ($FAILSTEP_RATE * 10)/1" | bc)
			done
		elif [ $i -eq 5 ]
	   	then
			FAILSTEP_RATE=$(echo "scale=3; 25/1" | bc)
		else
			FAILSTEP_RATE=$(echo "scale=3; 50/1" | bc)
		fi
		# Run through different k_off_pseudo rates: 0.5, 1, 5, 10, 50, 100
		for j in `seq 1 6`;
		do
			K_OFF_PSEUDO=$(echo "scale=2; 0.5/1" | bc)
			if [ $(($j % 2)) -eq 1 ]; then
				for pow in `seq 2 2 $j`;
				do
					K_OFF_PSEUDO=$(echo "scale=2; $K_OFF_PSEUDO * 10" | bc)
				done
			else
				K_OFF_PSEUDO=$(echo "scale=2; $K_OFF_PSEUDO * 2" | bc)
				for pow in `seq 3 2 $j`;
				do
					K_OFF_PSEUDO=$(echo "scale=2; $K_OFF_PSEUDO * 10" | bc)
				done
			fi
			# Run through desired MT lengths (2um, 4um, 6um, 8, 10um)
			for k in `seq 1 5`;
			do 
				MT_LENGTH=$(echo "scale=1; $k * 250" | bc)
				FILE_NAME="$K_OFF_RATIO"
				FILE_NAME+="_"
				if [ $i -le 3 ]; then
					FILE_NAME+="0"
				fi
				FILE_NAME+="$FAILSTEP_RATE"
				FILE_NAME+="_"
				if [ $j -le 1 ]; then
					FILE_NAME+="0"
				fi
				FILE_NAME+="$K_OFF_PSEUDO"
				FILE_NAME+="_"
				echo $FILE_NAME
				FILE_NAME+="$MT_LENGTH"
				echo RUNNING NEW SIM: filename is $FILE_NAME
				echo '(format is k-off-ratio _ failstep-rate _ k-off-pseudo _ sites)'
				# Update YAML files
				yq w -i params.yaml k_off_ratio $K_OFF_RATIO
				yq w -i params.yaml failstep_rate $FAILSTEP_RATE
				yq w -i params.yaml k_off_pseudo $K_OFF_PSEUDO
				yq w -i params.yaml length_of_microtubule $MT_LENGTH
				# Run sim for these parameter values
				./sim params.yaml $FILE_NAME
				# Remove unnecessary files after each sim (save mem.)
				MOTOR_FILE="$FILE_NAME"
				MOTOR_FILE+="_motorID.file"
				rm $MOTOR_FILE
				XLINK_FILE="$FILE_NAME"
				XLINK_FILE+="_xlinkID.file"
				rm $XLINK_FILE
				MT_FILE="$FILE_NAME"
				MT_FILE+="_MTcoord.file"
				rm $MT_FILE
			done
		done
	done
done	
echo END PARAMETER SCAN
