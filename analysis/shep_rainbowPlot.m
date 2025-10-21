
clear variables;

sim_name = 'out_final_newCombos/shep_0.75nM_30nM_8_1000_0.6kT_3x_5x_0';

start_frame = 1;
frames_per_plot = 200; % in n_datapoints; number of timesteps per output plot
end_frame = -1;  % set to -1 to run until end of data

% Load parameter structure
file_dir = '..';  % Default; only change if you move CyLaKS output files
params = load_parameters(sprintf('%s/%s', file_dir, sim_name));

% Open occupancy data file 
occupancy_filename = sprintf('%s/%s_occupancy.file', file_dir, sim_name);
occupancy = zeros(params.max_sites, params.n_mts, params.n_datapoints) - 1;
occupancy = load_data(occupancy, occupancy_filename, '*int');

xlink_speciesID = 1;
motor_speciesID = 2;

xlink_raw_data = occupancy; 
motor_raw_data = occupancy; 

xlink_raw_data(xlink_raw_data ~= xlink_speciesID) = 0;
xlink_raw_data(xlink_raw_data == xlink_speciesID) = 1;
motor_raw_data(motor_raw_data ~= motor_speciesID) = 0;
motor_raw_data(motor_raw_data == motor_speciesID) = 1;

xlink_avg_occupancy = zeros([params.max_sites params.n_mts]);
motor_avg_occupancy = zeros([params.max_sites params.n_mts]);

motor_avg_occupancy_tot = zeros([params.max_sites 1]);
xlink_avg_occupancy_tot = zeros([params.max_sites 1]);


% Read in and average occupancy data over all datapoints
i_plot = 1;
for i = 1:1:int32(params.n_datapoints)
    for i_pf = 1 : 1 : params.n_mts
        motor_avg_occupancy(:, i_pf) = motor_avg_occupancy(:, i_pf) + double(motor_raw_data(:, i_pf, i)) ./ frames_per_plot;
        xlink_avg_occupancy(:, i_pf) = xlink_avg_occupancy(:, i_pf) + double(xlink_raw_data(:, i_pf, i)) ./ frames_per_plot;
        motor_avg_occupancy_tot(:, 1) = motor_avg_occupancy_tot(:, 1) + double(motor_raw_data(:, i_pf, i)) ./ (frames_per_plot * params.n_mts);
        xlink_avg_occupancy_tot(:, 1) = xlink_avg_occupancy_tot(:, 1) + double(xlink_raw_data(:, i_pf, i)) ./ (frames_per_plot * params.n_mts);
    end
    if (mod(i, frames_per_plot) == 0)
        smooth_window = 32; % should be equivalent to diffraction limit
        motor_occupancy = smoothdata(motor_avg_occupancy, 'movmean', smooth_window);
        xlink_occupancy = smoothdata(xlink_avg_occupancy, 'movmean', smooth_window);
        motor_occupancy_tot = smoothdata(motor_avg_occupancy_tot, 'movmean', smooth_window);
        xlink_occupancy_tot = smoothdata(xlink_avg_occupancy_tot, 'movmean', smooth_window);  
        xlink_occu_vs_t(:, i_plot) = xlink_occupancy_tot; 
        i_plot = i_plot + 1;

        
        % Reset arrays to zero before we start counting again 
        for i_pf = 1 : 1 : params.n_mts
            motor_avg_occupancy(:, i_pf) = 0;
            xlink_avg_occupancy(:, i_pf) = 0;
        end
        motor_avg_occupancy_tot(:, 1) = 0;
        xlink_avg_occupancy_tot(:, 1) = 0;
    end
end
%}

xlink_occu_vs_t = flip(xlink_occu_vs_t, 1);

N = size(xlink_occu_vs_t, 2);
t_delta = params.time_per_datapoint * frames_per_plot;
colors = hsv(N);
distance = linspace(0, params.max_sites * params.site_size, params.max_sites);
figure()

plot3(distance, t_delta*ones(size(distance)), xlink_occu_vs_t(:, 1), 'Color', colors(1, :), 'LineWidth' ,2)
hold on
for i = 2 : N
    plot3(distance, i*t_delta*ones(size(distance)), xlink_occu_vs_t(:, i), 'Color', colors(i, :), 'LineWidth', 2)
end
grid on
xlim([0 8]);
xlabel('Distance (um)');
ylabel('Time (s)');
zlabel('MAP Occupancy');
set(gca, 'FontSize', 14);


