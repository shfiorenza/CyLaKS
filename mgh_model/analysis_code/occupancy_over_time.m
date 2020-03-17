clear all
% Dynamic variables
sim_name = 'test';
mt_lengths = [1000, 500];
% Pseudo-constant variables
n_steps = 40000000;
n_datapoints = 10000;
delta_t = 0.000025; 
species_ids = 1; % [1, 2];
species_labels = ["xlink", "motor"];
% Calculated quantities
n_mts = length(mt_lengths);
n_sites_max = max(mt_lengths);
n_species = length(species_ids);
time_per_datapoint = n_steps * delta_t / n_datapoints; 
% File stuff
file_directory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
occupancy_file_struct = '%s_occupancy.file';

occupancy_file_name = sprintf(file_directory, sprintf(occupancy_file_struct, sim_name));
occupancy_data_file = fopen(occupancy_file_name);
occupancy_data = fread(occupancy_data_file, n_mts * n_sites_max * n_datapoints, '*int');
fclose(occupancy_data_file);
occupancy_data = reshape(occupancy_data, n_sites_max, n_mts, n_datapoints);

fractional_occupancy = zeros(n_mts, n_datapoints, n_species);

for i_data = 1 : n_datapoints
    for i_mt = 1 : n_mts
        n_sites = mt_lengths(i_mt);
        for i_site = 1 : n_sites
            species_id = occupancy_data(i_site, i_mt, i_data);
            if species_id == 0
                continue;
            end
            fractional_occupancy(i_mt, i_data, species_id) = ...
                fractional_occupancy(i_mt, i_data, species_id) + 1.0 / double(n_sites);
        end
    end
end

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
