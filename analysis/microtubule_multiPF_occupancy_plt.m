clear variables;

sim_name = 'shep_multiPF_0_0.131_8';
sim_name = 'test2'

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
plot(xlink_avg_occupancy_tot, 'LineWidth', 2);

ylim([0 1]);


%{
for i_mt = 1 : n_mts
   for i_species = 1 : n_species
      fractional_occupancy(i_mt, :, i_species) = smoothdata(fractional_occupancy(i_mt, :, i_species)); 
   end
end

fig1 = figure();
hold all;

for i_mt = 1 : n_mts
    plot(linspace(0, delta_t * n_steps, n_datapoints), squeeze(fractional_occupancy(i_mt, :, :)));
end
legend(["static mt", "dynamic mt"]);
%}
