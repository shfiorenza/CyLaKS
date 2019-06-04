#!/usr/bin/python
import os
import numpy as np
import scipy
import logging
import subprocess 
import matlab.engine
from math import sqrt
from scipy.optimize import least_squares
from subprocess import call
from subprocess import Popen
MATLAB = matlab.engine.start_matlab()
np.set_printoptions(suppress=True)

sim_name = "least_squares"
param_name = "params_" + sim_name
log_file = sim_name + "_long_largeSteps" + ".scan"
mt_lengths = [2, 4, 6, 8, 10, 14]
exp_endtags = [1.1, 1.3, 1.45, 1.65, 1.8, 2.15]
exp_err_x = [0.8, 0.8, 0.8, 0.8, 0.8, 0.8]
exp_err_y = [0.1, 0.1, 0.1, 0.15, 0.15, 0.2]
exp_sigma = np.array(map(lambda x,y: sqrt(x**2 + y**2), exp_err_x, exp_err_y))
params_0 = np.array([8.75, 718, 117, 0.107])
labels = ["k_off_i", "k_off_i_stalled", "k_hydrolyze", "k_hydrolyze_stalled"]
param_bounds = ([1, 10, 15, 0.01], [50, 1000, 150, 15])
step_size = [5, 50, 10, 5]
call_no = 0

log = logging.getLogger()
log.setLevel(logging.DEBUG)
format = logging.Formatter('%(asctime)s - %(message)s', '%Y-%m-%d %H:%M:%S')
terminal_log_handler = logging.StreamHandler()
terminal_log_handler.setFormatter(format)
log.addHandler(terminal_log_handler)
file_log_handler = logging.FileHandler(log_file)
file_log_handler.setFormatter(format)
log.addHandler(file_log_handler)

def endtag_lengths(params):
    call('make clean-output', shell=True)
    global call_no
    iteration_no = int(call_no / 4) + 1
    sub_no = int(call_no % 4) + 1
    log.info('Beginning of iteration {}.{}'.format(iteration_no, sub_no))
    log.info('Params: {}'.format(params))
    # Copy base param file and alter parameters as desired
    call('cp params_base_endtag.yaml ' + param_name + ".yaml", shell=True)
    yaml_edit = "yq w -i " + param_name + ".yaml" + " motors."
    for i in range(len(params)):
        call(yaml_edit + labels[i] + " " + repr(params[i]), shell=True)

    # Make params and execution commands for sims with different MT lengths
    exe_commands = []
    param_files = []
    for i_length in reversed(range(len(mt_lengths))):
        n_sites = mt_lengths[i_length] * 125
        param_file = param_name + "_" + repr(n_sites) + ".yaml"
        call("cp " + param_name + ".yaml " + param_file, shell=True)
        yaml_edit = "yq w -i " + param_file + " microtubules."
        call(yaml_edit + "length[0] " + repr(n_sites), shell=True)
        cmd = "./sim " + param_file + " " + sim_name + "_" + repr(n_sites)
        exe_commands.append(cmd)
        param_files.append(param_file)

    sims = [ Popen(sim_exe, shell=True) for sim_exe in exe_commands ] 
    
    for sim in sims: sim.wait()
    for file in param_files: call("rm " + file, shell=True)
    weighted_errors = []
    for i_length in range(len(mt_lengths)):
        n_sites = mt_lengths[i_length] * 125
        name = sim_name + "_" + repr(n_sites)
        sim_endtag = MATLAB.get_endtag_length(str(sim_name), float(n_sites))
        err = exp_endtags[i_length] - sim_endtag
        weighted_err = err / exp_sigma[i_length]
        weighted_errors.append(weighted_err)

    log.info('Weighted errors: {}'.format(weighted_errors))
    call_no+=1
    return weighted_errors

already_made = os.path.isfile('./sim')
if not already_made: call("make CFG=release sim", shell=True)
log.info('Start of gradient descent parameter optimization')
log.info('Initial parameters: {}'.format(params_0))
res = least_squares(endtag_lengths, params_0, bounds=param_bounds, \
        diff_step=step_size, verbose=2, xtol=None)
