#!/usr/bin/python
import os
import numpy as np
import scipy
import logging
import subprocess
import matlab.engine
from scipy.optimize import least_squares
from subprocess import call
from subprocess import Popen
from math import sqrt
MATLAB = matlab.engine.start_matlab()
np.set_printoptions(suppress=True)

sim_name_base = "kif4a_coop_opt_summit"
param_file_base = "params/kif4a.yaml"
log_file = sim_name_base + ".scan"

param_label = ["neighb_neighb_energy", "gaussian_range",
               "gaussian_amp_solo", "gaussian_ceiling_bulk"]
param_initialVal = [1.6, 1000, 0.95, 4.75]
param_bounds = ([0.5, 100, 0.1, 1.5], [3.0, 1500, 1.5, 6])
step_size = [0.1, 100, 0.1, 0.1]

# param_label = ["interaction_energy", "lattice_coop_range",
#               "lattice_coop_Emax_solo", "lattice_coop_Emax_bulk",
#               "sigma_off_ii", "sigma_ATP"]
#param_initialVal = [0.75, 850, 0.75, 3.75, 0.175, 2.3]
#param_bounds = ([0.5, 100, 0.1, 1.5, 0.05, 0.05], [3.0, 1500, 1.5, 4.5, 5, 5])
#step_size = [0.1, 100, 0.1, 0.1, 0.1, 0.5]

n_runs = [50, 125, 200, 300, 400, 500]

# Kif4A concentrations in pM
kif4a_conc = [20, 50, 80, 120, 220, 420]
# Experimental Kif4A run lengths at each concentration
exp_runlengths = [0.9735, 1.310, 2.420, 1.659, 1.964, 2.855]
exp_err_runlengths = [0.18, 0.32, 0.35, 0.94, 0.31, 0.72]
# Experiemntal Kif4A life times at each concentration
exp_lifetimes = [1.821, 2.083, 7.096, 5.233, 8.308, 17.95]
exp_err_lifetimes = [0.56, 0.75, 1.8, 5.9, 2.6, 3.9]
# Experimental Kif4A velocities at each concentration
exp_velocities = [598.6, 709.9, 360.8, 311.2, 308.52, 180.7]
exp_err_velocities = [76, 110, 49, 78, 40, 38]

# Create logger to record history of optimizer
log = logging.getLogger()
log.setLevel(logging.DEBUG)
format = logging.Formatter('%(asctime)s - %(message)s', '%Y-%m-%d %H:%M:%S')
# Add handler that prints log to terminal
terminal_log_handler = logging.StreamHandler()
terminal_log_handler.setFormatter(format)
log.addHandler(terminal_log_handler)
# Add handler that prints log to file
file_log_handler = logging.FileHandler(log_file)
file_log_handler.setFormatter(format)
log.addHandler(file_log_handler)

call_no = 0


def kif4a_coop_scaling(params):
    global call_no
    # Keep track of which iteration of the gradient descent algorithm we're on
    iteration_no = int(call_no / len(params))
    # To avoid overwriting data, keep track of sub-iterations too (one parameter is varied at a time)
    sub_no = int(call_no % len(params))
    log.info("Beginning of iteration {}.{}".format(iteration_no, sub_no))
    log.info("Parameters: {}".format(params))
    # Delete all output files that weren't saved (not enough mem to keep them all)
    call("./clean_output.sh", shell=True)
    # Make names, param files, and execution commands for sims with different kif4a concentrations
    sim_names = []
    param_files = []
    exe_commands = []
    # Set up order so that largest sim (highest conc) runs first -- not sure if this actually matters
    for i_conc in range(len(kif4a_conc)):
        conc = kif4a_conc[i_conc]
        # Generate unique simulation name for this iteration, sub-iteration, and kif4a concentration
        sim_name = sim_name_base + "_" + \
            repr(iteration_no) + "." + repr(sub_no) + "_" + repr(conc)
        sim_names.append(sim_name)
        # Name parameter file in the same fashion
        param_file = "temp_params_" + sim_name + ".yaml"
        # Copy base parameter file and make desired parameter edits
        call("cp " + param_file_base + " " + param_file, shell=True)
        yaml_edit = "yq eval -i "
        call(yaml_edit + '".motors.c_bulk = ' + repr(conc * 0.001) + '" ' + param_file, shell=True)  # pM to nM
        call(yaml_edit + '".motors.n_runs_to_exit = ' + repr(n_runs[i_conc]) + '" ' + param_file, shell=True)
        # FIXME 
        for i in range(len(params)):
            call(yaml_edit + '".motors.' + param_label[i] + ' = '
                 " " + repr(params[i]) + '" ' + param_file, shell=True)
        # Add parameter file to param_files array so we can rm them later
        param_files.append(param_file)
        # Generate command to execute this simulation
        # exe_cmd = "./cylaks.exe " + param_file + " " + sim_name
        exe_cmd = "singularity exec --bind $PWD cylaks_latest.sif cylaks.exe  " + param_file + " " + sim_name
        # Add command to exe_commands array so we can launch them all at once
        exe_commands.append(exe_cmd)
    # Launch all sims (with different kif4a concentrations) at once
    sims = [Popen(sim_exe, shell=True) for sim_exe in exe_commands]
    # Wait for all simulations to finish
    for sim in sims:
        sim.wait()
    # Calculate errors
    weighted_errors = []
    for i_conc in range(len(kif4a_conc)):
        kif4a_stats = MATLAB.get_motor_stats(str(sim_names[i_conc]))
        log.info("For sim {}:".format(sim_names[i_conc]))
        log.info("  Measured stats: {}".format(kif4a_stats))
        err_runlength = (
            exp_runlengths[i_conc] - kif4a_stats[0][0]) / exp_err_runlengths[i_conc]
        log.info("      Runlength error: {}".format(err_runlength))
        err_lifetime = (
            exp_lifetimes[i_conc] - kif4a_stats[0][2]) / exp_err_lifetimes[i_conc]
        log.info("      Lifetime error: {}".format(err_lifetime))
        err_velocity = (
            exp_velocities[i_conc] - kif4a_stats[0][4]) / exp_err_velocities[i_conc]
        log.info("      Velocity error: {}".format(err_velocity))
        # Only use lifetime and velocity so as to not overconstrain 
        weighted_errors.append(err_lifetime**2 + err_velocity**2)
    log.info("Weighted errors: {}".format(weighted_errors))
    # Move log file into output folder to keep a record of all runs
    call("mv *.log grad_descent_output", shell=True)
    # Remove temporary parameter files
    for file in param_files:
        call("rm " + file, shell=True)
    call_no += 1
    return weighted_errors


#already_made = os.path.isfile("./cylaks.exe")
#if not already_made:
#    call("./install.sh", shell=True)
ready_for_output = os.path.isfile("/grad_descent_output/")
if not ready_for_output:
    call("mkdir grad_descent_output", shell=True)
log.info("Start of gradient descent parameter optimization")
log.info("Initial parameters: {}".format(param_initialVal))
res = least_squares(kif4a_coop_scaling, param_initialVal,
                    bounds=param_bounds, diff_step=step_size,
                    x_scale='jac', verbose=2)  # , xtol=None)
