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

sim_name_base = "kif4a_coop_optimization"
param_file_base = "params_processivity.yaml"
log_file = sim_name_base + ".scan"

param_label = ["lattice_coop_amp", "lattice_coop_range"]
param_initialVal =  np.array([7, 500])
param_bounds = ([1, 1], [100, 1500])
step_size = [0.5, 10]

# Kif4A concentrations in pM
kif4a_conc = [20, 50, 80, 120, 220, 420]
# Experimental Kif4A run lengths at each concentration
exp_runlengths = [0.97, 1.31, 2.42, 1.66, 1.96, 2.86]
exp_err_runlengths = [0.18, 0.32, 0.35, 0.94, 0.31, 0.72]
# Experiemntal Kif4A life times at each concentration
exp_lifetimes = [1.8, 2.1, 7.1, 5.2, 8.3, 17.9]
exp_err_lifetimes = [0.6, 0.7, 1.7, 5.9, 2.6, 3.9]
# Experimental Kif4A velocities at each concentration
exp_velocities = [600, 710, 360, 310, 310, 180]
exp_err_velocities = [80, 110, 50, 80, 40, 40]

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
    # Make names, param files, and execution commands for sims with different kif4a concentrations
    sim_names = []
    param_files = []
    exe_commands = []
    for conc in kif4a_conc:
        # Generate unique simulation name for this iteration, sub-iteration, and kif4a concentration
        sim_name = sim_name_base + "_" + repr(iteration_no) + "." + repr(sub_no) + "_" + repr(conc)
        sim_names.append(sim_name)
        # Name parameter file in the same fashion
        param_file = "temp_params_" + sim_name + ".yaml"
        # Copy base parameter file and make desired parameter edits
        call("cp " + param_file_base + " " + param_file, shell=True)
        yaml_edit = "yq w -i " + param_file + " motors."
        for i in range(len(params)): 
            call(yaml_edit + param_label[i] + " " + repr(params[i]), shell=True)
        # Add parameter file to param_files array so we can rm them later
        param_files.append(param_file)
        # Generate command to execute this simulation
        exe_cmd = "./sim" + param_file + " " + sim_name
        # Add command to exe_commands array so we can launch them all at once 
        exe_commands.append(exe_cmd)
    # Launch all sims (with different kif4a concentrations) at once
    sims = [ Popen(sim_exe, shell=True) for sim_exe in exe_commands ] 
    # Wait for all simulations to finish
    for sim in sims: sim.wait()
    # Remove temporary parameter files
    for file in param_files: call("rm " + file, shell=True)
    # Move log file into output folder to keep a record of all runs 
    call("mv *.log grad_descent_output", shell=True)
    # Delete all other output files (not enough mem to keep them all)
    call("make clean-output", shell=True)
    # Calculate errors
    weighted_errors = []
    for i_conc in range(len(kif4a_conc)):
        kif4a_stats = MATLAB.get_motor_stats(str(sim_names[i]))
        err_runlength = (exp_runlengths[i_conc] - kif4a_stats[0]) / exp_err_runlengths[i_conc]
        err_lifetime = (exp_lifetimes[i_conc] - kif4a_stats[1]) / exp_err_lifetimes[i_conc]
        err_velocity = (exp_velocities[i_conc] - kif4a_stats[2]) / exp_err_velocities[i_conc]
        weighted_errors.append(err_runlength + err_lifetime + err_velocity)
    log.info('Weighted errors: {}'.format(weighted_errors))
    call_no+=1
    return weighted_errors

already_made = os.path.isfile('./sim')
if not already_made: call("make -j12 CFG=release sim", shell=True)
ready_for_output = os.path.isfile('/grad_descent_output/')
if not ready_for_output: call("mkdir grad_descent_output", shell=True)
log.info('Start of gradient descent parameter optimization')
log.info('Initial parameters: {}'.format(param_initialVal))
res = least_squares(kif4a_coop_scaling, param_initialVal, bounds=param_bounds,\
        diff_step=step_size, verbose=2, xtol=None)
