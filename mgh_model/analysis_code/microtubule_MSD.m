clear all;
simName = "MT_diffu_test";
n_sites = [1000, 500];
max_sites = max(n_sites);
n_mts = length(n_sites);
starting_point = 001;
max_tau = 1.5;  % in seconds
% Pseudo-constant variables
delta_t = 0.000025;
n_steps = 6000000;
n_datapoints = 10000;
time_per_datapoint = delta_t * n_steps / n_datapoints;
active_datapoints = n_datapoints - starting_point;
site_size = 0.008; % in um

fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
mtFileStruct = '%s_mt_coord.file';
mtFileName = sprintf(fileDirectory, sprintf(mtFileStruct, simName));
mt_data_file = fopen(mtFileName);
mt_raw_data = fread(mt_data_file, n_mts * n_datapoints, '*double');
fclose(mt_data_file);
mt_data = reshape(mt_raw_data, n_mts, n_datapoints);

starting_tau = time_per_datapoint;
tau_increment = time_per_datapoint;
n_taus = max_tau / tau_increment;

n_entries = zeros(n_mts, n_taus);
CSD = zeros(n_mts, n_taus);
MSD = zeros(n_mts, n_taus);

for i_tau=1:n_taus  
    tau = starting_tau + (i_tau - 1) * tau_increment;
    tau_step = int32(tau / time_per_datapoint);
    for i_mt=1:n_mts
        for i_data = (tau_step + starting_point):tau_step:n_datapoints
            prev_pos = mt_data(i_mt, i_data - tau_step);
            cur_pos = mt_data(i_mt, i_data);
            dist = (cur_pos - prev_pos) * site_size;
            dist_sq = dist * dist;
            CSD(i_mt, i_tau) = CSD(i_mt, i_tau) + dist_sq; 
            n_entries(i_mt, i_tau) = n_entries(i_mt, i_tau) + 1;
        end
        MSD(i_mt, i_tau) = CSD(i_mt, i_tau) / n_entries(i_mt, i_tau);  
    end
end

D = zeros(n_mts, 1);

for i_mt=1:n_mts
    % Use basic linear regression to find slope
    y = MSD(i_mt, :)';
    x = zeros(length(MSD(i_mt, :)), 1);
    for i_tau = 1:1:n_taus
        x(i_tau) = starting_tau + (i_tau - 1) * tau_increment;
    end
    X = [ones(length(x), 1) x];
    m = X\y;
    D(i_mt) = m(2) / 2
end

plot(linspace(starting_tau, max_tau, n_taus), MSD);