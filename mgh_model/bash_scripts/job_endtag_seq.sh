echo "Running on $(hostname --fqdn)"

#for K_OFF_RATIO in {10}; do
K_OFF_RATIO=40;
	K_OFF_I_ST=$(echo "scale=2; 7.2 * $K_OFF_RATIO" | bc)
#	for K_HYDRO in {90}; do
	K_HYDRO=105;
		#for JAM_RATIO in {900}; do
		JAM_RATIO=900;
        	K_HYDRO_ST=$(echo "scale=3; $K_HYDRO / $JAM_RATIO" | bc)
#			mkdir jam_ratio_$JAM_RATIO
			for MT_LENGTH in {250,500,750,1000,1250,1750}; do
				PARAM_FILE="params_endtag_"
				PARAM_FILE+=$K_OFF_RATIO
				PARAM_FILE+="_"
				PARAM_FILE+=$K_HYDRO
				PARAM_FILE+="_"
				PARAM_FILE+=$JAM_RATIO
				PARAM_FILE+="_"
				PARAM_FILE+=$MT_LENGTH
				PARAM_FILE+=".yaml"
				echo PARAM FILE is $PARAM_FILE
				cp params_base_endtag.yaml $PARAM_FILE
                yq w -i $PARAM_FILE motors.k_off_i_stalled $K_OFF_I_ST
                yq w -i $PARAM_FILE motors.k_hydrolyze $K_HYDRO
                yq w -i $PARAM_FILE motors.k_hydrolyze_stalled $K_HYDRO_ST
                yq w -i $PARAM_FILE microtubules.length[0] $MT_LENGTH
				SIM_NAME="endtag_"
				SIM_NAME+=$K_OFF_RATIO
				SIM_NAME+="_"
				SIM_NAME+=$K_HYDRO
				SIM_NAME+="_"
				SIM_NAME+=$JAM_RATIO
				SIM_NAME+="_"
				SIM_NAME+=$MT_LENGTH
				SIM_NAME+="_long"
#				SIM_NAME+="_medlow_kOn"
				echo SIM NAME is $SIM_NAME
				./sim $PARAM_FILE $SIM_NAME
				echo Deleted $PARAM_FILE
				rm $PARAM_FILE
				echo Moved $SIM_NAME to jam_ratio_$JAM_RATIO
#				mv $SIM_NAME* jam_ratio_$JAM_RATIO
			done
#		done
#	done
#done
