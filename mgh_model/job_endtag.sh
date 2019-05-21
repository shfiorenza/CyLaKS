#!/bin/bash

# === List of SBATCH arguments ===
#SBATCH --account=ucb-summit-smr	# allocations account
#SBATCH --nodes=4				# number of cluster nodes, abbreviated by -N
#SBATCH --ntasks=96				# number of parallel process
#SBATCH --time=15:00:00			# walltime, abbreviated by -t
#SBATCH --qos=normal        	# quality of service/queue
#SBATCH --job-name=et_scan			# name of job
#SBATCH --output=et_scan-%j.out    # name of the stdout redirection file
#SBATCH --error=et_scan-%j.err		# name of the stderr redirection file

# === Purge all modules and load needed ones ===
module purge
module load intel
module load gsl
module list
echo "Current variables CC=${CC} and CXX=${CXX}"

# === Run the program ===
echo "Running on $(hostname --fqdn)"

#for K_HYDRO in {75,90,105,120,135}; do
K_HYDRO=60;
	for JAM_RATIO in {100,200,300,400,500,600,700,800}; do
		for MT_LENGTH in {250,500,750,1000,1250,1750}; do
			PARAM_FILE="params_endtag_"
			PARAM_FILE+=$K_HYDRO
			PARAM_FILE+="_"
			PARAM_FILE+=$JAM_RATIO
			PARAM_FILE+="_"
			PARAM_FILE+=$MT_LENGTH
			PARAM_FILE+=".yaml"
			echo PARAM FILE is $PARAM_FILE
			SIM_NAME="endtag_"
			SIM_NAME+=$K_HYDRO
			SIM_NAME+="_"
			SIM_NAME+=$JAM_RATIO
			SIM_NAME+="_"
			SIM_NAME+=$MT_LENGTH
			echo SIM NAME is $SIM_NAME
#			srun -N 1 -n 1 ./sim $PARAM_FILE $SIM_NAME &
			./sim $PARAM_FILE $SIM_NAME
		done
	done
#done
wait

mkdir scan_output
#for K_HYDRO in {75,90,105,120,135}; do
K_HYDRO=60;
	mkdir k_hydrolyze_$K_HYDRO
	for JAM_RATIO in {100,200,300,400,500,600,700,800}; do
		mkdir jam_ratio_$JAM_RATIO
		for MT_LENGTH in {250,500,750,1000,1250,1750}; do
			SIM_NAME="endtag_"
			SIM_NAME+=$K_HYDRO
			SIM_NAME+="_"
			SIM_NAME+=$JAM_RATIO
			SIM_NAME+="_"
			SIM_NAME+=$MT_LENGTH
			echo moved $SIM_NAME to jam_ratio_$JAM_RATIO
			mv $SIM_NAME* jam_ratio_$JAM_RATIO
		done
		echo moved jam_ratio_$JAM_RATIO to k_hydrolyze_$K_HYDRO
		mv jam_ratio_$JAM_RATIO k_hydrolyze_$K_HYDRO
	done
	echo moved k_hydrolyze_$K_HYDRO to scan_output
	mv k_hydrolyze_$K_HYDRO scan_output
#done
