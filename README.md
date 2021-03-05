# CyLaKS
The **Cy**toskeleton **La**ttice **K**intetic **S**imulator
## Installation
To begin, simply clone this repo:
```
git clone https://github.com/Betterton-Lab/CyLaKS
cd CyLaKS
```
### Dependencies 
CyLaKS requires the following libraries to run:
	* [gsl] (http://www.gnu.org/software/gsl/)
	* [yaml-cpp] (https://github.com/jbeder/yaml-cpp)

On most Ubuntu distributions, you can install these dependences as follows:
```	
sudo apt-get install libgsl-dev
sudo apt-get install libyaml-cpp-dev
```

In order to run any of the included scripts, [yq](https://github.com/mikefarah/yq) is also required. There are a variety of ways to install this utility, so consult the documentation to find the most convenient means for your local environment. 
### Building from source 
To use the provided installation script, you must have CMake (version 3.14+) installed. On most Ubuntu distributions, CMake can be installed via:
```
sudo apt-get install cmake
```

Once CMake is installed, simply call the installation script:
```
./install.sh 
```

There are additional flags to build CyLaKS in debug mode or clean the build directory. To view these options, add the '-h' flag:
```
./install.sh -h
```
Once compiled, see CyLaKS run options by calling the executable:
```
./sim 
```
## Running CyLaKS
Input file types 
 .yaml -> contains parameters; can be created/editing using any text-editor
Output file types
 .log  -> contains a record of simulation outputs; access using any text-editor
 .file -> contains raw data from sims; use included MATLAB code to analyze
### Demos
Different demos can be found in the demos/ folder. 
### Custom simulations
placeholder