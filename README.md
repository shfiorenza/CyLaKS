File suffixes correspond to the following:
 .yaml -> contains parameters; can be created/editing using any text-editor
 .log  -> contains a record of sim outputs; access using any text-editor
 .file -> contains raw data from sims; use included MATLAB code to analyze

Included makefile has various options:

# RELEASE mode (optimized compiler flags; MUCH faster runtime)
	make -jN  					(N = number of cores available to your computer)
# DEBUG mode
    make -jN CFG=debug			(N = number of cores available to your computer)
# Once compiled, see run options by running
	./sim 
# Compiling on summit (make sure to use COMPILE node)
	module purge
	module load intel gsl
	make -j12 LOC=summit 
# To run on Summit, use one of the included job files:
	sbatch job_xyz.sh
# To check job status, use squeue command: 
	squeue --user=YOU
