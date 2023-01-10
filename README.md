# CyLaKS
The **Cy**toskeleton **La**ttice-based **K**inetic **S**imulator
## Installation
### Using pre-built Docker image
A pre-built image of CyLaKS is available as a [Docker](https://www.docker.com/) image. To download this image, run
```
docker pull shfiorenza/cylaks
```
Depending your local Docker permissions, you may need to add `sudo` before the command. If this is your first time using Docker, 
you may also have to initialize the Docker Daemon. On most Linux distributions, you can simply run
```
dockerd
```
On Windows and MacOS, you can use [Docker Desktop](https://www.docker.com/products/docker-desktop). 

Once you have downloaded the CyLaKS image, you can use the provided `launch_docker.sh` script to easily initialize the Docker environemnt and run simulations within it. (Data files are still be output to your local directory.) To see a full list of script run options, use the _help_ flag: 
```
./launch_docker.sh -h
```
Before running any simulations, the Docker container must first be appropriately launched in the background. To do this, use the _initialize_ flag:  
```
./launch_docker.sh -i
```
Again, you may need to add 'sudo' depending on your local Docker permissions. Once the container is initialized, you can start the interactive simulation launcher using the _run_ flag:
```
./launch_docker.sh -r
```
You can also manually launch CyLaKS simulations via Docker while the container is running:
```
docker exec -it cylaks cylaks.exe [parameter-file] [sim-name] [optional-flags]
```
See the _Running Simulations_ section for more information on the interactive launcher and quick-launch syntax. 

To close the cylaks container, use the `stop` command:
```
docker stop cylaks
```
This will leave the docker image on your computer. To fully remove the docker image, use the _clean up_ flag: 
```
./launch_docker.sh -c
```
Finally, if the Dockerhub image does not work for whatever reason, try building the environment and simulation images locally on your machine:
```
./launch_docker.sh -e
./launch_docker.sh -b
```

### Building from source 
If you plan on compiling or instally CyLaKS locally, first clone the repo:
```
git clone https://github.com/Betterton-Lab/CyLaKS
cd CyLaKS
```
In order to properly compile, CyLaKS requires the following libraries:
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
There are additional flags to build CyLaKS in debug mode or clean the build directory. To view these options, use the _help_ flag:
```
./install.sh -h
```
Once compiled, you can run the interactive CyLaKS launcher by simply calling the executable:
```
./cylaks.exe
```
## Running simulations
With no provided input arguments, the 'cylaks.exe' executable will run an interactive simulation launcher. 
To launch simulations more quickly, use the _quick-launch syntax_ with the `cylaks.exe` executable as follows:
```
./cylaks.exe [parameter-file] [sim-name] [optional-flags]
```
The parameter file must be a YAML file. See the files included in `params/` for an example of proper formatting. 

The sim name sets the prefix of all output files and can be whatever you'd like. 

The optional flags allow you to choose from and configure a variety of test/demo modes, as discussed below. 
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
