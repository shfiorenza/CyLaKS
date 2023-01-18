# CyLaKS
The **Cy**toskeleton **La**ttice-based **K**inetic **S**imulator
## Installation
### Using a pre-built Docker image
A pre-built image of CyLaKS is available as a [Docker](https://www.docker.com/) image. To download this image, run
```
docker pull shfiorenza/cylaks
```
Depending on your local Docker permissions, you may need to add `sudo` before the command. If this is your first time using Docker, 
you may also have to initialize the Docker Daemon. On most Linux distributions, you can run
```
dockerd
```
On Windows and MacOS, you can use [Docker Desktop](https://www.docker.com/products/docker-desktop). 

Once you have downloaded the CyLaKS image, you can use the provided `launch_docker.sh` script to easily initialize the Docker environment and run simulations within it. (Data files are still output to your local directory.) 

To see a full list of script run options, use the _help_ (-h) flag: 
```
./launch_docker.sh -h
```
Before running simulations, the Docker container must first be appropriately launched in the background. To do this, use the _initialize_ (-i) flag:  
```
./launch_docker.sh -i
```
Again, you may need to add 'sudo' depending on your local Docker permissions. Once the container is initialized, you can start the interactive simulation launcher using the _run_ (-r) flag:
```
./launch_docker.sh -r
```
You can also manually launch CyLaKS simulations via Docker while the container is running:
```
docker exec -it cylaks cylaks.exe [parameter-file] [sim-name] [optional-flags]
```
See the _Running simulations_ section for more information on the interactive launcher and launch syntax. 

To close the cylaks container, use the `stop` command:
```
docker stop cylaks
```
This will leave the Docker image on your computer for later use. 

To fully remove the Docker image, use the _clean up_ (-c) flag: 
```
./launch_docker.sh -c
```
Finally, if the Dockerhub image does not work for whatever reason, try building the images locally:
```
./launch_docker.sh -e
./launch_docker.sh -b
```

### Building from source 
If you plan on compiling or installing CyLaKS locally, first clone the repo:
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
Once CMake is installed, call the installation script to compile CyLaKS. An executable named `cylaks.exe.` will be added to the main folder. 
```
./install.sh 
```
There are additional flags to build CyLaKS in debug mode or clean the build directory. 

To view these options, use the _help_ (-h) flag:
```
./install.sh -h
```
Once compiled, you can run the interactive CyLaKS launcher by simply calling the executable:
```
./cylaks.exe
```
See the _Running simulations_ section for more information. 
## Running simulations
With no other provided input arguments, the 'cylaks.exe' executable will run an interactive simulation launcher. 

If you choose to run a test/demo mode, the simulation name and parameters will be automatically generated. 
You may be prompted to input additional parameters that are specific to the chosen mode. 

If you choose to run a 'regular' simulation, you simply need to input the desired parameter file and simulation name. 
The parameter file must be a YAML file and have values set for all system parameters, even if not all are used. 
See the files included in `params/` for an example of proper formatting. 
The simulation name sets the prefix of all data output files and can be whatever you'd like. 

To launch simulations more quickly, you can use the following _quick-launch_ syntax:
```
./cylaks.exe [parameter-file] [sim-name]
```
You can also manually run test/demo modes with your own simulation name by using additional input arguments
```
./cylaks.exe [parameter-file] [sim-name] [test mode] [test mode args]
```
To see more information on available demos and test modes, see the _Demos_ and _Test modes_ sections. 
### Demos
#### Motility of mutant heterodimeric motors
Mutant heterodimers are modeled as motors with two distinct heads: one 'normal' head that steps unidirectionally, and one 'mutant' head that diffuses randomly in addition to regular stepping. As the diffusion rate of the second head is increased, the heterodimeric motor becomes less efficient at stepping as a whole. 

You can choose this demo by selecting the `kinesin_mutant` mode in the interactive launcher. Alternatively, you can launch it manually:
```
./cylaks.exe [params] [sim-name] kinesin_mutant
```
This mode does not accept any additonal arguments. To change the diffusion of the second head, modify the 'd_i' variable in the 'xlinks' section of the parameter file. 
#### Motors stepping along a heterogeneous tubulin lattice
Each binding site is a distinct object within CyLaKS, allowing us to easily modify the properties of each. In this demo mode, a certain fraction of tubulin binding sites are modified to have lower binding affinity relative to normal sites. The so-called 'mutant' sites are randomly distributed throughout the tubulin lattice, and can lead to interesting traffic effects when motors are introduced. 

You can choose this demo by selecting the 'hetero_tubulin' mode in the interactive launcher. Alternatively, you can launch it manually:
```
./cylaks.exe [params] [sim-name] hetero_tubulin [p_mutant] [bind_aff]
```
where `p_mutant` is the probability of each site to be a mutant (ranges from 0 to 1) and `bind_aff` is the binding affinity of a regular site relative to mutant sites. 
E.g., setting `bind_aff = 3` will make motors bind to regular sites 3x as often and unbind 3x less often. 
#### Kif4A end-tag formation and response to filament ablation
#### Microtubule sliding assays to probe force response of crosslinkers
These modes allow the user to directly insert a certain number of doubly bound crosslinkers into anti-parallel microtubule overlaps, circumventing the need to wait for the system to equilibrate over time. 

The first mode fixes the x-position of each microtubule and examines the average separation between them as a function of crosslinker number. You can choose this demo by selecting the 'filament_separation' mode in the interactive launcher. Alternatively, you can launch it manually: 
```
./cylaks.exe [params] [sim-name] filament_separation [n_xlinks]
```
where `n_xlinks` is the number of doubly bound crosslinkers to insert at the beginning of the simulation. 

The second mode examines the force response of crosslinkers to two different types of externally imposed forces. One type is meant to reflect optical trapping experiments, in which a constant velocity is imposed on one microtubule while the other microtubule is held fixed. The second type is meant to reflect kinesin gliding experiments, in which both microtubules have a variable gliding velocity and move in opposite directions. You can choose this demo by selecting the 'filament_forced_slide' mode in the interactive launcher. Alternatively you can launch it manually: 
```
./cylaks.exe [params] [sim-name] filament_forced_slide [n_xlinks] [slide_vel] [force_flag] 
```
or
```
./cylaks.exe [params] [sim-name] filament_forced_slide [n_xlinks] [slide_vel] [force_flag] [t_pause] [pause_dur]
```
where `n_xlinks` is the number of doubly bound crosslinkers to insert at the beginning of the simulation, `slide_vel` is the imposed velocity of the microtubule(s) in nm/s, and 'force_flag' selects between the two different imposed force types. Setting `force_flag = 0` selects the optical trapping mode, and the input velocity will be the constant velocity at which one microtubule is dragged. Setting `force_flag = 1` selects the kinesin gliding mode, and the input velocity will be the constant gliding velocity _in the absence of any other force_ (i.e, crosslinkers can slow them down). The two optional inputs `t_pause` and `pause_dir` (in seconds) allow for a temporary pause in force application (usually used in force mode 0) to see how the system responds after being allowed to mechanically reorganize after induced sliding.

### Parameters
#### General
```
seed                    Seed of random number generator
kbT                     Boltzmann constant multiplied by temperature; pN*nm
eta                     Viscosity of liquid; (pN/um^2)*s
dt                      Time elapsed each kmc-bd step; s
t_run                   Time to run post-equilibration; s
t_equil                 Minimum time to equilibrate; s
t_snapshot              Time between each data output; s
dynamic_equil_window    Set to <=0 to disable dynamic equil; s
verbosity               Verbosity level; 0 (quiet) to 3 (max)    
```
#### Filaments
```
count                   Number of filaments in simulation
n_subfilaments          For multiple protofilaments in a MT
periodic_barrel         Whether or not barrel wraps around completelty
radius                  Radius of rod (or barrel for MTs); nm
site_size               Length of each binding site; nm
n_bd_per_kmc            BD iterations done per KMC iteration
n_sites                 Length of each filament; n_sites
polarity                0 (1) sets plus-end to i=0 (n_sites - 1)
x_initial               Starting x-coord of filament COM; nm
y_initial               Starting y-coord of filament COM; nm
x_immobile_until        Time when fil. can move in x; s
y_immobile_until        Time when fil. can move in y; s
rotation_enabled        Toggle rotation within  x-y plane
wca_potential_enabled   Toggle WCA potential between filaments
f_applied               Force applied to filament COM; pN
```
#### Motors
```
n_runs_to_exit          Number of post-equil. runs to trigger an exit
gaussian_range          Range of Gaussian interaction; n_sites
gaussian_amp_solo       Amplitude of interaction for one motor; kBT
gaussian_ceiling_bulk   Ceiling of amp. for many motors; kBT
gaussian_stepping_coop  Do long-rane FX apply to stepping?
neighb_neighb_energy    Short-range interaction magnitude; kBT
t_active                Time at which motors/ATP flows in; s
k_on                    Binding rate of ADP heads to MT; 1/(nM*s)
c_bulk                  Bulk concentration of motors; nM
c_eff_bind              For 2nd ADP head when 1st is bound; nM
k_on_ATP                ATP binding rate; 1/(micromolar*s)
c_ATP                   Bulk concentration of ATP; micromolar
k_hydrolyze             Rate that ATP->ADPP occurs in motors; 1/s
k_off_i                 Unbinding rate for ADPP-bound heads; 1/s
k_off_ii                Unbinding rate for ADPP-bound heads; 1/s
applied_force           Perpetually applied force on motors; pN
internal_force          Internal necklinker tension; pN
sigma_ATP               ATP-binding distance parameter; nm
sigma_off_i             Singly-bound off-rate distance parameter; nm
sigma_off_ii            Doubly-bound off-rate distance parameter; nm
endpausing_active       Toggles if motor step off plus-ends
tethers_active          Toggles tethers (Kif4A stalks or 'tails')
k_tether                Tethering rate; 1/(nM*s)
c_eff_tether            Effective concentration of free_teth motors
k_untether              Untethering rate when bound; 1/s
r_0                     Rest length of stalk (or tether); nm
k_spring                Spring constant of tether; pN/nm
k_slack                 '' but when shorter than rest length
```
#### Crosslinkers 
```
neighb_neighb_energy    Energy between two PRC1 neighbors
t_active                Time at which crosslinkers flow in; s
k_on                    Bulk binding rate;  1/(nM*s)
c_bulk                  Bulk concentration; nM
c_eff_bind              Effective conc. of 2nd head binding; nM
k_off_i                 Unbinding rate while singly-bound;  1/s
k_off_ii                Unbinding rate while doubly-bound;  1/s
d_i                     Diffusion coefficient; um^2/s
d_ii                    Diffusion coefficient; um^2/s
d_side                  Diffusion coefficient; um^2/s (side-stepping)
r_0                     Rest length of coiled-coil domain; nm
k_spring                Spring constant of CC-domain; pN/nm
theta_0                 Rotational rest angle
k_rot                   Rotational spring constant;
p_diffuse_off_end       Chance to diffuse off plus or minus ends of MTs
```

## Analyzing output
CyLaKS has two types of output files:
 * .log: plain-text record of simulation output
 * .file: raw data collected during simulation
Log files can be viewed in any text editor, e.g., notepad. Data files must be parsed using the MATLAB scripts found in `analysis/` 

To use the included analysis scripts, simply change the `sim_name` variable to the name of the simulation you wish to analyze. Note that this should be the raw sim name without directory or file suffixes, e.g., identical to the name of the log file without the `.log` suffix. 

## Codebase development
### Adding new parameters 
### Test modes
#### Doubly bound crosslinker diffusion
#### Doubly bound crosslinker binding
#### Long-range binding interaction between motors
#### Long-range stepping interaction between motors
