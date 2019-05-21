File suffixes correspond to the following:
 .yaml -> contains parameters; can be created/editing using any text-editor
 .log  -> contains a record of sim outputs; access using any text-editor
 .file -> contains raw data from sims; use included MATLAB code to analyze

Included makefile has various options:
# DEBUG mode
    make sim
# RELEASE mode (optimized compiler flags; MUCH faster runtime)
	make CFG=release sim
# Compiling on summit (make sure to use COMPILE node)
	module purge
	module load intel gsl
	make CFG=release LOC=summit sim
# To run on Summit, use one of the included job files:
	sbatch job_xyz.sh
# To check job status, use squeue command: 
	squeue --user=YOU
