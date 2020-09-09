# name of output folder for all scan runs
SCAN_NAME="kif4a_coop"
# name of slurm file that will launch each run
JOB_SCRIPT="job_gradOpt_lattice.slurm"
# desired directory for output folder
DIR="/scratch/summit/shfi4480/overlap_analysis/mgh_model/"

SCAN_OUT="${DIR}${SCAN_NAME}/"
# if output folder doesn't exist, create it
if [ ! -d "${SCAN_OUT}" ]; then
	echo "Creating scan output folder @ ${SCAN_OUT}"
	mkdir ${SCAN_OUT}
fi

PARAM_ONE="ProteinConfig.yaml proteins[0]freeLength";
#PARAM_TWO="RunConfig.yaml rngSeed";

for VAL_ONE in 100 # 7 8
do
#	for VAL_TWO in 1 2 3 4 5
#	do
	#	for VAL_THREE in 0.01 0.1 1
	#	do	
			# name of output folder for this particular scan run
	#		NAME="r_${VAL_ONE}_${VAL_TWO}_${VAL_THREE}"
	#		NAME="r_${VAL_ONE}_${VAL_TWO}"
			NAME="r_${VAL_ONE}"
			echo " ***** INITIALIZING ${NAME} ***** "
			OUTPUT="${SCAN_OUT}${NAME}"
			# if specific output already exists, do not overwrite it in any circumstance
			if [ -d "${OUTPUT}" ]; then
				echo "    ** already exists; skipping **"
			# otherwise, copy necessary files and modify parameters/job_name
			else
				mkdir ${OUTPUT}
				# copy executable and other auxiliary files
				python3 scripts/CopyToRun.py ${OUTPUT}
				# copy parameter files
				rsync -avP RunConfig.yaml ${OUTPUT}
				rsync -avP ProteinConfig.yaml ${OUTPUT}
				# if initialization files exist, copy them too
				if [ -f "TubuleInitial.dat" ]; then
					rsync -avP TubuleInitial.dat ${OUTPUT}
				fi
				if [ -f "ProteinInitial.dat" ]; then
					rsync -avP ProteinInitial.dat ${OUTPUT}
				fi
				# copy job script 
				rsync -avP ${JOB_SCRIPT} ${OUTPUT}
				# navigate to newly-created directory
				cd ${OUTPUT}
				# edit params
				yq w -i ${PARAM_ONE} ${VAL_ONE}e-3
			#	yq w -i ${PARAM_TWO} $((1234+VAL_TWO))
			#	yq w -i ${PARAM_THREE} [${VAL_THREE} ${VAL_THREE}]
				# edit job name
				REGEX="s/--job-name=.*/--job-name=${NAME}/g;"
				perl -i -p -e ${REGEX} ${JOB_SCRIPT}
				# launch job
				sbatch ${JOB_SCRIPT}
				# return to previous directory where script was called
				cd - 
			fi
	#	done
#	done
done
