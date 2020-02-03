clear all;
simName = "test";
n_mts = 2;
max_tau = 300; 
% Pseudo-constant variables
site_size = 0.008; % in um
delta_t = 2.5e-3;
n_steps = 6000000;
n_datapoints = 10000;
time_per_datapoint = delta_t * n_steps / n_datapoints;

fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
mtFileStruct = '%s_mt_coord.file';
mtFileName = sprintf(fileDirectory, sprintf(mtFileStruct, simName));
mt_data_file = fopen(mtFileName);
mt_raw_data = fread(mt_data_file, n_mts * n_datapoints, '*int');
fclose(mt_data_file);
mt_data = reshape(mt_raw_data, n_mts, n_datapoints);

min_tau = time_per_datapoint;
tau_increment = time_per_datapoint;
n_taus = (max_tau - min_tau) / tau_increment;
max_tau_step = max_tau / time_per_datapoint;

n_entries = zeros(n_mts, n_taus);
CSD = zeros(n_mts, n_taus);
MSD = zeros(n_mts, n_taus);

for i_tau=1:n_taus  
    tau = min_tau + (i_tau - 1) * tau_increment;
    tau_step = tau / time_per_datapoint;
    for i_mt=1:1:n_mts
        for i_data = max_tau_step + 1:n_datapoints
            prev_pos = mt_data(i_mt, i_data - int32(tau_step));
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
y_int = zeros(n_mts, 1);

for i_mt=1:n_mts
    % Use basic linear regression to find slope
    y = MSD(i_mt, :)';
    x = zeros(length(MSD(i_mt, :)), 1);
    for i_tau = 1:1:n_taus
        x(i_tau) = min_tau + (i_tau - 1) * tau_increment;
    end
    X = [ones(length(x), 1) x];
    m = X\y;
    D(i_mt) = m(2) / 2;
    y_int(i_mt) = m(1);
    fprintf("MT %i: D = %g, y_int = %g\n", i_mt, D(i_mt), y_int(i_mt));
end
figure();
plot(linspace(min_tau, max_tau, n_taus), MSD);
%figure();
%plot(linspace(0, n_datapoints, n_datapoints), double(mt_data)*site_size);