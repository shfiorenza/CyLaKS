%{
clear all;
% Dynamic variables
base_name = 'prc1_coop';
mt_lengths = [500];
coop_energies = [2.25]; 
xlink_concs = [0.9056, 2:2:40]; 
%xlink_concs = [0.9056, 38];
seeds = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
% Experimental data 
exp_scaling = [1.0, 1.39, 10.2, 36.6, 78.1, 168.3];
exp_scaling_err = [0.05, 0.05, 0.2230, 0.1039, 0.1216, 0.1016]; % fractional
exp_concs = [0.9056, 1.802, 9.251, 18.49, 27.74, 37.03];
%{
exp_scaling = [0.09, 0.14, 0.68, 2.3, 5.0, 11.00];
exp_scaling = exp_scaling/min(exp_scaling);
exp_concs = [0.92, 1.8, 9.2, 19, 28, 38];
%}
% Pseudo-constant variables
n_datapoints = 10000;
xlink_speciesID = 1;
n_mts = length(mt_lengths);
n_sites_max = max(mt_lengths);
n_energies = length(coop_energies);
n_concs = length(xlink_concs);
n_seeds = length(seeds);
% File stuff
fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
fileStructure = '%s_occupancy.file';

norm_intensity = zeros([n_energies n_concs n_seeds n_datapoints]);
for i_energy = 1 : n_energies
    for i_conc = 1 : n_concs
        for i_seed = 1 : n_seeds
            sim_name = base_name;
            if n_energies > 1
                sim_name = sprintf("%s_g%", sim_name, coop_energies(i_energy));
            end
            if n_concs > 1
                if i_conc ~= 1
                    sim_name = sprintf("%s_%i", sim_name, xlink_concs(i_conc));
                else
                    sim_name = sprintf("%s_%i", sim_name, 0);
                end
            end
            if n_seeds > 1
                sim_name = sprintf("%s_%i", sim_name, seeds(i_seed));
            end
            fileName = sprintf(fileDirectory, sprintf(fileStructure, sim_name));
            xlink_data_file = fopen(fileName);
            xlink_data = fread(xlink_data_file, n_mts * n_sites_max * n_datapoints, '*int');
            fclose(xlink_data_file);
            xlink_data = reshape(xlink_data, n_sites_max, n_mts, n_datapoints);
            
            xlink_data(xlink_data ~= xlink_speciesID) = 0;
            xlink_data(xlink_data == xlink_speciesID) = 1;
            
            for i_mt = 1 : n_mts
                n_sites = mt_lengths(i_mt);
                for i_data = 1 : n_datapoints
                    occupancy = xlink_data(:, i_mt, i_data);
                    mt_intensity = sum(occupancy);
                    norm_intensity(i_energy, i_conc, i_seed, i_data) = mt_intensity / n_sites;
                end
            end
        end
    end
end

norm_intensity_avg = zeros(n_energies, n_concs);
norm_intensity_err = zeros(n_energies, n_concs);

for i_energy = 1 :n_energies
    for i_conc = 1 : n_concs
        seed_avgs = zeros(n_seeds, 1);
        for i_seed = 1 : n_seeds
            seed_avgs(i_seed) = mean(norm_intensity(i_energy, i_conc, i_seed, 2500:10000));
        end
        norm_intensity_avg(i_energy, i_conc) = mean(seed_avgs);
        variance = 0.0;
        for i_seed = 1 : n_seeds
           seed_avg =  mean(norm_intensity(i_energy, i_conc, i_seed, 2500:10000));
           diff_sq = (norm_intensity_avg(i_energy, i_conc) - seed_avg)^2;
           variance = variance + diff_sq / double(n_seeds - 1);
        end
        norm_intensity_err(i_energy, i_conc) = sqrt(variance / n_seeds);
    end
end
%}
fig1 = figure();
hold on
avg_I0 = mean(norm_intensity_avg(:, 1));
exp_prediction = avg_I0 * exp_scaling;
exp_prediction_err = exp_prediction.*exp_scaling_err;
for i_energy = 1 : n_energies
    errorbar(xlink_concs, norm_intensity_avg(i_energy, :), norm_intensity_err(i_energy, :), ...
        '.', 'LineWidth', 2, 'MarkerSize', 14);
end
errorbar(exp_concs, exp_prediction, exp_prediction_err, '--s', 'LineWidth', 1.5, 'MarkerSize', 7);

%title('Interaction energy of -2.5 kBT for PRC1-PRC1 neighbors');
xlabel('PRC1 Concentration (nM)', 'FontSize', 14);
ylabel('Normalized Intensity (A.U.)', 'FontSize', 14);
legend_label = ["Simulation", "Experiment"];
%{
energy_labels = strcat('E_{int} = -',strtrim(cellstr(num2str(coop_energies','%.2f'))), ' k_BT');
for i=length(legend_label):length(coop_energies)
    legend_label(i + 1) = energy_labels(i);
end
%}
legend(legend_label, 'location', 'northwest', 'FontSize', 14);