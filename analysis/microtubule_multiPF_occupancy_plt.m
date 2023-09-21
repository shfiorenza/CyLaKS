clear variables;

sim_name = 'shep_10nM_20nM_8_0.0kT';
sim_name = 'coop_131_3.6nM_0nM_8_1.45kT_0';
sim_name = 'endtags_3/endtag_0.0524_25_0.1nM_20nM_8_1.375kT_500_3';

% Load parameter structure
file_dir = '../shepherding_baseline';  % Default; only change if you move CyLaKS output files
file_dir = '../coop_fit3_10xDiff/usable_output';
file_dir = '..';
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

xlink_avg_occupancy = zeros([params.n_datapoints params.n_mts]);
motor_avg_occupancy = zeros([params.n_datapoints params.n_mts]);

motor_avg_occupancy_tot = zeros([params.n_datapoints 1]);
xlink_avg_occupancy_tot = zeros([params.n_datapoints 1]);

for i_data = 1 : params.n_datapoints
    for i_mt = 1 : params.n_mts
        n_sites = params.mt_lengths(i_mt);
        for i_site = 1 : n_sites
            species_id = occupancy(i_site, i_mt, i_data);
            if species_id == 0
                continue;
            elseif species_id == xlink_speciesID
                xlink_avg_occupancy(i_data, i_mt) = xlink_avg_occupancy(i_data, i_mt) + 1.0 / double(n_sites);
                xlink_avg_occupancy_tot(i_data, 1) = xlink_avg_occupancy_tot(i_data, 1) + 1.0 / double(n_sites * params.n_mts);
            elseif species_id == motor_speciesID
                motor_avg_occupancy(i_data, i_mt) = motor_avg_occupancy(i_data, i_mt) + 1.0 / double(n_sites);
            end
        end
    end
end

fig = figure();

plot(xlink_avg_occupancy);
hold on
plot(motor_avg_occupancy); 
plot(xlink_avg_occupancy_tot, 'LineWidth', 2);

ylim([0 1]);

saveas(fig, 'yeet.jpg', 'jpg');

