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

sim_base = "least_squares_scan"
params_base = "params_endtag.yaml"
log_file = sim_base + "_teth_singleMT" + ".scan"
mt_lengths = [2, 4, 6, 8, 10, 14]
exp_endtags_0 = [1.1, 1.3, 1.5, 1.7, 1.8, 2.2]
exp_endtags_1 = [1.3, 1.6, 1.8, 2.0, 2.3, 2.8]
exp_endtags_4 = [1.6, 2.1, 2.7, 3.2, 3.8, 4.9]
exp_err_x = [0.8, 0.8, 0.8, 0.8, 0.8, 0.8]
exp_err_y_0 = [0.1, 0.1, 0.1, 0.15, 0.15, 0.2]
exp_err_y_1 = [0.1, 0.6, 0.2, 0.30, 0.30, 0.2]
exp_err_y_4 = [0.1, 0.1, 0.1, 0.15, 0.15, 0.2]
exp_sigma_0 = np.array(map(lambda x,y: sqrt(x**2 + y**2), exp_err_x, exp_err_y_0))
exp_sigma_1 = np.array(map(lambda x,y: sqrt(x**2 + y**2), exp_err_x, exp_err_y_1))
exp_sigma_4 = np.array(map(lambda x,y: sqrt(x**2 + y**2), exp_err_x, exp_err_y_4))

initial_params =  np.array([0.1, 0.05, 5000])
labels = ["k_tether", "k_untether", "c_eff_tether"]
param_bounds = ([0.00001, 0.00005, 50], [15, 15, 15000])
step_size = [0.01, 0.01, 100]
xlink_concs = [1, 4]
#params_0 = np.array([8.75, 718, 117, 0.107])
#labels = ["k_off_i", "k_off_i_stalled", "k_hydrolyze", "k_hydrolyze_stalled"]
#param_bounds = ([1, 10, 15, 0.01], [50, 1000, 150, 15])
#step_size = [5, 50, 10, 5]

log = logging.getLogger()
log.setLevel(logging.DEBUG)
format = logging.Formatter('%(asctime)s - %(message)s', '%Y-%m-%d %H:%M:%S')
terminal_log_handler = logging.StreamHandler()
terminal_log_handler.setFormatter(format)
log.addHandler(terminal_log_handler)
file_log_handler = logging.FileHandler(log_file)
file_log_handler.setFormatter(format)
log.addHandler(file_log_handler)

call_no = 0
def endtag_lengths(params):
    call('make clean-output', shell=True)
    global call_no
    iteration_no = int(call_no / len(params))
    sub_no = int(call_no % len(params))
    log.info('Beginning of iteration {}.{}'.format(iteration_no, sub_no))
    log.info('Params: {}'.format(params))
    # Make params and execution commands for sims with different MT lengths
    exe_commands = []
    param_files = []
    # Copy base param file and alter parameters as desired
    for conc in xlink_concs:
        sim_name = sim_base + "_" + repr(iteration_no) + "." + repr(sub_no) + "_" + repr(conc)
        param_file = "params_" + sim_base + "_" + repr(conc) + ".yaml"
        call("cp " + params_base + " " + param_file, shell=True)
        call("yq w -i " + param_file + " xlinks.c_bulk " + repr(float(conc)/10), shell=True)
        yaml_edit = "yq w -i " + param_file + " motors."
        for i in range(len(params)): call(yaml_edit + labels[i] + " " + repr(params[i]), shell=True)
        for i_length in reversed(range(len(mt_lengths))):
            n_sites = mt_lengths[i_length] * 125
            temp_params = "params_temp_" + repr(conc) + "_" + repr(n_sites) + ".yaml"
            call("cp " + param_file + " " + temp_params, shell=True)
            yaml_edit = "yq w -i " + temp_params + " microtubules."
            call(yaml_edit + "length[0] " + repr(n_sites), shell=True)
            cmd = "./sim " + temp_params + " " + sim_name + "_" + repr(n_sites)
            exe_commands.append(cmd)
            param_files.append(temp_params)

        call("rm " + param_file, shell=True)

    sims = [ Popen(sim_exe, shell=True) for sim_exe in exe_commands ] 
    for sim in sims: sim.wait()
    for file in param_files: call("rm " + file, shell=True)

    weighted_errors = []
    for conc in xlink_concs:
        sim_name = sim_base + "_" + repr(iteration_no) + "." + repr(sub_no) + "_" + repr(conc)
        for i_length in range(len(mt_lengths)):
            n_sites = mt_lengths[i_length] * 125
            sim_endtag = MATLAB.get_endtag_length(str(sim_name), float(n_sites))
            err = 0;
            weighted_err = 0;
            if conc == 1:
                err = exp_endtags_1[i_length] - sim_endtag
                weighted_err = err / exp_sigma_1[i_length]
            elif conc == 4:
                err = exp_endtags_4[i_length] - sim_endtag
                weighted_err = err / exp_sigma_4[i_length]
            else:
                call("WHAT")
            weighted_errors.append(weighted_err)

    call("mv *.log grad_descent_output", shell=True);
    log.info('Weighted errors: {}'.format(weighted_errors))
    call_no+=1
    return weighted_errors

already_made = os.path.isfile('./sim')
if not already_made: call("make CFG=release sim", shell=True)
ready_for_output = os.path.isfile('/grad_descent_output/')
if not ready_for_output: call("mkdir grad_descent_output", shell=True);
log.info('Start of gradient descent parameter optimization')
log.info('Initial parameters: {}'.format(initial_params))
res = least_squares(endtag_lengths, initial_params, bounds=param_bounds,\
        diff_step=step_size, verbose=2, xtol=None)
