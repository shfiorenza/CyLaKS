File suffixes correspond to the following:
 .yaml -> contains parameters; can be created/editing using any text-editor
 .file -> contains raw data from sims; use included MATLAB code to analyze
 .log  -> contains a record of sim outputs; access using any text-editor

For makefile:
# DEBUG mode
    make sim
# RELEASE mode (optimized compiler flags; MUCH faster runtime)
	make CFG=release sim
# Compiling on summit; set LOC=summit (Make sure to use COMPILE node)
	module purge
	module load intel gsl
	make CFG=release LOC=summit sim
# To run on summit, use one of the included job files:
	sbatch job_xyz.sh
