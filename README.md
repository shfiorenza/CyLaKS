# CyLaKS

The **Cy**toskeleton **La**ttice **K**intetic **S**imulator

# Installation

Simply clone the git repo:
```
git clone https://github.com/Betterton-Lab/CyLaKS
cd CyLaKS
```

# Dependencies 

CyLaKS requires the following libraries to run:
	* [yaml-cpp] (https://github.com/jbeder/yaml-cpp)
	* [gsl] (http://www.gnu.org/software/gsl/)

On most Ubuntu distributions, you can install these dependences as follows:
```
	sudo apt-get install libyaml-cpp-dev
	sudo apt-get install libgsl-dev
```

In order to run any of the included scripts, [yq](https://github.com/mikefarah/yq) is also required. There are a variety of ways to install this utility, so consult the documentation to find the most convenient means for your local environment. 



# Building from source 

Included makefile has various options:
# RELEASE mode (optimized compiler flags; MUCH faster runtime)
	make -jN  					(N = number of cores available to your computer)
# DEBUG mode
    make -jN CFG=debug			(N = number of cores available to your computer)
# Demos

Different demos can be found in the demos/ folder. 

# Running new simulations

Input file types 
 .yaml -> contains parameters; can be created/editing using any text-editor
Output file types
 .log  -> contains a record of simulation outputs; access using any text-editor
 .file -> contains raw data from sims; use included MATLAB code to analyze


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
