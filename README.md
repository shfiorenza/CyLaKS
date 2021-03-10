# CyLaKS
The **Cy**toskeleton **La**ttice-based **K**inetic **S**imulator
## Installation
To begin, clone the repo:
```
git clone https://github.com/Betterton-Lab/CyLaKS
cd CyLaKS
```
### Dependencies 
CyLaKS requires the following libraries to run:
 * [gsl](http://www.gnu.org/software/gsl/)
 * [yaml-cpp](https://github.com/jbeder/yaml-cpp)

On most Ubuntu distributions, you can install these dependences as follows:
```	
sudo apt-get install libgsl-dev
sudo apt-get install libyaml-cpp-dev
```
In order to run most of the included bash scripts found in `scripts/`, [yq](https://github.com/mikefarah/yq) is also required. There are a variety of ways to install this utility, so consult the documentation to find the most convenient means for your local environment. 
### Building from source 
To use the provided installation script, you must have CMake (version 3.14+) installed. On most Ubuntu distributions, CMake can be installed via:
```
sudo apt-get install cmake
```

Once CMake is installed, call the installation script to compile CyLaKS. An executable will be added to the main folder. 
```
./install.sh 
```

There are additional flags to build CyLaKS in debug mode or clean the build directory. To view these options, add the '-h' flag:
```
./install.sh -h
```
Once compiled, you can observe CyLaKS run options by calling the executable:
```
./cylaks
```
## Running simulations
To run CyLaKS, use the `cylaks` executable as follows:
```
./cylaks [parameter-file] [sim-name] [optional-flags]
```
The parameter file must be a YAML file. See the files included in `params/` for an example of proper formatting. 

The sim name sets the prefix of all output files and can be whatever you'd like. 

### Demos
To select from currently available demos, use the `run_demos.sh` script. 
### Test modes
placeholder

## Analyzing output
CyLaKS has two types of output files:
 * .log: plain-text record of simulation output
 * .file: raw data collected during simulation
Log files can be viewed in any text editor, e.g., notepad. Data files must be parsed using the MATLAB scripts found in `analysis/` 

To use the included analysis scripts, simply change the `sim_name` variable to the name of the simulation you wish to analyze. Note that this should be the raw sim name without directory or file suffixes, e.g., identical to the name of the log file without the `.log` suffix. 
