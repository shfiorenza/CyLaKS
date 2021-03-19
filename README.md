# CyLaKS
The **Cy**toskeleton **La**ttice-based **K**inetic **S**imulator
## Installation
To begin, clone the repo:
```
git clone https://github.com/Betterton-Lab/CyLaKS
cd CyLaKS
```
### Using pre-built Docker image
A pre-built image of CyLaKS is available as a [Docker](https://www.docker.com/) image. 

To download this image, run
```
docker pull shfiorenza/cylaks
```
Depending your local Docker permissions, you may need to add `sudo` to the command above. 

If this is the first time using Docker on your machine, you may have to initialize the Docker Daemon. On Linux distributions, run
```
dockerd
```
On Windows and MacOS, you can use [Docker Desktop](https://www.docker.com/products/docker-desktop). 

Once you have the CyLaKS image, you can launch a Docker container named `cylaks` in the background via the provided script:
```
./launch_docker.sh
```
Again, you may need to add 'sudo' depending on your local Docker permissions. 

Once the container is running, you can launch CyLaKS simulations via Docker:
```
docker exec -it cylaks cylaks.exe [parameter-file] [sim-name] [optional-flags]
```
You can also select from available demos by using the 
### Building from source 
CyLaKS requires the following libraries:
 * [gsl](http://www.gnu.org/software/gsl/)
 * [yaml-cpp](https://github.com/jbeder/yaml-cpp)

On most Linux distributions, you can install these dependences as follows:
```	
sudo apt-get install libgsl-dev
sudo apt-get install libyaml-cpp-dev
```
In order to use most of the included bash scripts found in `scripts/`, [yq](https://github.com/mikefarah/yq) is also required. There are a variety of ways to install this utility, so consult the documentation to find the most convenient means for your local environment. 

To use the provided installation script, you must have CMake (version 3.13+) installed. On most Linux distributions, CMake can be installed via:
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
./cylaks.exe
```
## Running simulations
To run CyLaKS, use the `cylaks.exe` executable as follows:
```
./cylaks.exe [parameter-file] [sim-name] [optional-flags]
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
