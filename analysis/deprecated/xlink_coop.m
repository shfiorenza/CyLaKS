
clear all;
% Dynamic variables
n_datapoints = 10000;
n_equil = 5000; 
base_name = 'prc1_coop_%snM_8_%skT_1.3x_%i'; 
fileDirectory = '/home/shane/projects/CyLaKS/out_coop8/%s';
mt_lengths = [2500, 2500, 2500, 2500, 2500, 2500, 2500, 2500];
coop_energies = [1.15]; % , 1.18];
e_labels = ["1.15"]; %, "1.18"]; 
%coop_energies = [0.9, 1.0, 1.1, 1.2, 1.3];
%e_labels = ["0.9", "1.0", "1.1", "1.2", "1.3"];
xlink_concs = [0.906, 1.8, 3.6, 7.2, 9.2, 13.6, 18.5, 27.7, 37.0];
c_labels = ["0.906", "1.8", "3.6", "7.2", "9.2", "13.6", "18.5", "27.7", "37.0"];
seeds = [0,1,2,3];
% Experimental data
exp_scaling = [1.0, 1.39, 10.2, 36.6, 78.1, 168.3];
exp_scaling_err = [0.05, 0.05, 0.2230, 0.1039, 0.1216, 0.1016]; % fractional
exp_concs = [0.9056, 1.802, 9.251, 18.49, 27.74, 37.03];
%exp_concs = exp_concs / 10.0;
sim_concs = xlink_concs; %[0.09, 0.18, 0.93, 1.8, 2.8, 3.7];
%{
exp_scaling = [0.09, 0.14, 0.68, 2.3, 5.0, 11.00];
exp_scaling = exp_scaling / min(exp_scaling);
exp_concs = [0.92, 1.8, 9.2, 19, 28, 38];
%}
% Pseudo-constant variables

xlink_speciesID = 1;
n_mts = length(mt_lengths);
n_sites_max = max(mt_lengths);
n_energies = length(coop_energies);
n_concs = length(xlink_concs);
n_seeds = length(seeds);
% File stuff
fileStructure = '%s_occupancy.file';

norm_intensity = zeros([n_energies n_concs n_seeds n_datapoints]);

for i_energy = 1:n_energies
    for i_conc = 1:n_concs
        for i_seed = 1:n_seeds
            sim_name = sprintf(base_name, c_labels(i_conc), e_labels(i_energy), seeds(i_seed))
           
            fileName = sprintf(fileDirectory, sprintf(fileStructure, sim_name));
            xlink_data_file = fopen(fileName);
            xlink_data = fread(xlink_data_file, n_mts * n_sites_max * n_datapoints, '*int');
            fclose(xlink_data_file);
            xlink_data = reshape(xlink_data, n_sites_max, n_mts, n_datapoints);

            xlink_data(xlink_data ~= xlink_speciesID) = 0;
            xlink_data(xlink_data == xlink_speciesID) = 1;

            for i_mt = 2:n_mts-1
                n_sites = mt_lengths(i_mt);

                for i_data = 1:n_datapoints
                    occupancy = xlink_data(:, i_mt, i_data);
                    mt_intensity = sum(occupancy);
                    mt_intensity_norm = mt_intensity / n_sites; 
                    norm_intensity(i_energy, i_conc, i_seed, i_data) = norm_intensity(i_energy, i_conc, i_seed, i_data) ...
                        + mt_intensity_norm / (n_mts-2);
                end

            end

        end

    end

end

norm_intensity_avg = zeros(n_energies, n_concs);
norm_intensity_err = zeros(n_energies, n_concs);

for i_energy = 1:n_energies
    for i_conc = 1:n_concs
        seed_avgs = zeros(n_seeds, 1);
        for i_seed = 1:n_seeds
            seed_avgs(i_seed) = mean(norm_intensity(i_energy, i_conc, i_seed, n_equil:n_datapoints));
        end

        norm_intensity_avg(i_energy, i_conc) = mean(seed_avgs);
        variance = 0.0;

        for i_seed = 1:n_seeds
            seed_avg = mean(norm_intensity(i_energy, i_conc, i_seed, n_equil:n_datapoints));
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
exp_prediction_err = exp_prediction .* exp_scaling_err;

errorbar(exp_concs, exp_prediction, exp_prediction_err, '--s', 'LineWidth', 2, 'MarkerSize', 10);

for i_energy = 1:n_energies
    errorbar(xlink_concs, norm_intensity_avg(i_energy, :), norm_intensity_err(i_energy, :), ...
        '.', 'LineWidth', 2, 'MarkerSize', 20);
end

%}
%title('Interaction energy of -2.5 kBT for PRC1-PRC1 neighbors');
xlabel('PRC1 Concentration (nM)', 'FontSize', 18);
ylabel('Normalized Intensity (A.U.)', 'FontSize', 18);
set(gca, 'FontSize', 18);
legend_label = ["Exp prediction", compose('Sim - %g kT', e_labels)];
%{
energy_labels = strcat('E_{int} =- ',strtrim(cellstr(num2str(coop_energies', '%.2f'))), ' k_BT');

for i = length(legend_label):length(coop_energies)
    legend_label(i + 1) = energy_labels(i);
end

%}
legend(legend_label, 'location', 'northwest', 'FontSize', 18);
