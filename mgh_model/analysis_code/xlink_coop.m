
clear variables;
baseName = 'coop_bind_longB';
mt_length = [1000];
E_int = [2.25]; %, 2.50, 2.75, 3.00, 3.25];
exp_scaling = [0.09, 0.14, 0.68, 2.3, 5.0, 11.00];
exp_scaling = exp_scaling/min(exp_scaling);
xlink_concs = [0.9, 1.8, 9.2, 19, 28, 38];
%xlink_concs = [1,2,10,20,28,38];

n_datapoints = 10000;
starting_point = 0;
xlink_speciesID = 1;
active_datapoints = n_datapoints - starting_point;
n_mts = length(mt_length);
n_sites_max = max(mt_length);
n_energies = length(E_int);
n_concs = length(xlink_concs);
norm_intensity = zeros([n_energies n_concs n_datapoints]);

fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
fileStructure = '%s_occupancy.file';

for i_energy=1:n_energies
    for i_conc=1:n_concs
        simName = sprintf('%s_%i_%i', baseName, int32(E_int(i_energy)*100), int32(xlink_concs(i_conc)*10))
        fileName = sprintf(fileDirectory, sprintf(fileStructure, simName));
        xlink_data_file = fopen(fileName);
        xlink_raw_data = fread(xlink_data_file, n_mts * n_sites_max * n_datapoints, '*int');
        fclose(xlink_data_file);
        xlink_data = reshape(xlink_raw_data, n_sites_max, n_mts, n_datapoints);
        
        xlink_data(xlink_data ~= xlink_speciesID) = 0;
        xlink_data(xlink_data == xlink_speciesID) = 1;
        
        for i_data=1:n_datapoints
            for i_mt=1:n_mts
                occupancy = xlink_data(:, i_mt, i_data);
                n_sites = mt_length(i_mt);
                for i_site=1:1:n_sites
                    norm_intensity(i_energy, i_conc, i_data) = ...
                        norm_intensity(i_energy, i_conc, i_data) + (double(occupancy(i_site)) / n_sites);
                end
            end
        end
    end
end

%norm_intensity = smoothdata(norm_intensity);
intensity_avgs = zeros([n_energies n_concs]);
for i_energy=1:n_energies
    for i_conc=1:1:n_concs
        intensity_avgs(i_energy,i_conc) =  mean(norm_intensity(i_energy, i_conc, 5000:10000));
    end
end

fig1 = figure();
hold on
avg_I0 = mean(intensity_avgs(:, 1));
exp_prediction = avg_I0 * exp_scaling;
plot(xlink_concs, exp_prediction, '^', 'LineWidth', 2, 'MarkerSize', 14);

for i_energy=1:n_energies
    plot(xlink_concs, intensity_avgs(i_energy, :), '*', 'LineWidth', 2, 'MarkerSize', 14);
end

%title('Interaction energy of -2.5 kBT for PRC1-PRC1 neighbors');
xlabel('PRC1 Concentration (nM)', 'FontSize', 14);
ylabel('Intensity (A.U.)', 'FontSize', 14);
legend_label = {'Experiment'}; %, 'No cooperativity'};
energy_labels = strcat('E_{int} = -',strtrim(cellstr(num2str(E_int','%.2f'))), ' k_BT');
for i=length(legend_label):length(E_int)
    legend_label(i + 1) = energy_labels(i);
end
legend(legend_label, 'location', 'northwest', 'FontSize', 14);