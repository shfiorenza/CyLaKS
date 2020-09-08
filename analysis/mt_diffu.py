import numpy as np
import matplotlib.pyplot as plt

n_mts = 4
max_tau = 30
simName = "test"
delta_t = 2.5e-5
n_steps = 6000000
n_datapoints = 10000
time_per_datapoint = n_steps * delta_t / n_datapoints
site_size = 0.008

fileDirectory = "/home/shane/Projects/overlap_analysis/mgh_model/%s"
mtFileStruct = "%s_mt_coord.file"
mtFilePath = fileDirectory % mtFileStruct % simName
mt_data = np.fromfile(mtFilePath, dtype=np.int32).reshape(n_datapoints, n_mts)

min_tau = time_per_datapoint
tau_increment = time_per_datapoint
tau_values = np.arange(min_tau, max_tau, tau_increment)

max_tau_step = int(np.amax(tau_values) / time_per_datapoint)
n_samples = n_datapoints - max_tau_step

CSD = np.zeros((n_mts, len(tau_values)))
for i_tau, tau in enumerate(tau_values):
    tau_step = int(tau / time_per_datapoint)
    for i_mt in range(n_mts):
        for i_data in np.arange(tau_step, n_datapoints):
            prev_pos = mt_data[i_data - tau_step, i_mt]
            cur_pos = mt_data[i_data, i_mt]
            dist = (prev_pos - cur_pos) * site_size
            dist_sq = dist**2
            CSD[i_mt, i_tau] += dist_sq

MSD = CSD / n_samples
diffusion_coeffs = np.zeros((n_mts, 1))
y_intercept = np.zeros((n_mts, 1))

A = np.vstack([tau_values, np.ones(len(tau_values))]).T
for i_mt in range(n_mts):
    m, c = np.linalg.lstsq(A, MSD[i_mt], rcond=None)[0]
    diffusion_coeffs[i_mt] = m / 2
    y_intercept[i_mt] = c
    print("For mt %i: D = %g, c = %g" % (i_mt, m/2, c))

plt.plot(tau_values, MSD[0])
plt.plot(tau_values, MSD[1])
plt.show()
# print(diffusion_coeffs)
