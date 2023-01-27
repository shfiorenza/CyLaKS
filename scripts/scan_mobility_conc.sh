#!/bin/bash
SCAN_NAME="motor_mobility"
echo Starting ${SCAN_NAME} scan
PARAM_FILE="params/kif4a.yaml"
echo Base parameter file is ${PARAM_FILE}

MOT_CONC=(0.02 0.05 0.08 0.120 0.220 0.420)
CONC_SCALE=(20 50 80 120 220 420)
N_RUNS=(50 125 200 300 500 800)
BASE_SEED=198261346419
MIN_SEED=0
MAX_SEED=3

N_THREADS=12
I_THREAD=0

echo "Need to run a batch of 24 sims. How many CPU threads to use?" 
read N_THREADS

for I_CONC in 0 1 2 3 4 5
do
	for (( SEED_NO=MIN_SEED; SEED_NO <= MAX_SEED; SEED_NO++))
	do
		FILE_NAME="${SCAN_NAME}_${CONC_SCALE[I_CONC]}_${SEED_NO}"
		TEMP_PARAMS="params_temp_${FILE_NAME}.yaml"
		cp $PARAM_FILE $TEMP_PARAMS
        yq eval -i ".seed = $(( $BASE_SEED + $SEED_NO ))" ${TEMP_PARAMS}
		# yq eval -i ".t_equil = ${T_EQUIL[I_CONC]}" ${TEMP_PARAMS}
        yq eval -i ".motors.c_bulk = ${MOT_CONC[I_CONC]}" ${TEMP_PARAMS}			
        yq eval -i ".motors.n_runs_to_exit = ${N_RUNS[I_CONC]}" ${TEMP_PARAMS}
		# Run sim for these parameter values
		echo Running new sim: file name is ${FILE_NAME}
#		./sim $TEMP_PARAMS $FILE_NAME &
		./cylaks.exe $TEMP_PARAMS $FILE_NAME &
		# singularity exec --bind $PWD cylaks_latest.sif cylaks.exe $TEMP_PARAMS $FILE_NAME &
		((I_THREAD++))
		if [ $I_THREAD -ge $N_THREADS ]
		then
			wait
			((I_THREAD=0))
		fi
	done
done
wait
rm params_temp_${SCAN_NAME}*
