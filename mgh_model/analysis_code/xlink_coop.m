clear variables;
%simName = 'coop1';
mt_length = [5000];
n_mts = length(mt_length);
n_sites_max = max(mt_length);
n_steps = 100000000;
delta_t = 0.00001;
n_datapoints = 10000;
starting_point = 0;
active_datapoints = n_datapoints - starting_point;
xlink_speciesID = 1;

fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
fileStructure = '%s_occupancy.file';

%exp_scaling = [1, 1.22, 6.65, 23.6, 50.5, 109.2];
exp_scaling = [0.0916, 0.112, 0.609, 2.15, 4.62, 10.00];
exp_scaling = exp_scaling/min(exp_scaling);
%exp_scaling = [1, 5.44, 19.26, 41.24, 89.25];
xlink_concs = [1.1, 2.1, 10.6, 21.8, 32.5, 43.1];
%xlink_concs = [2.1, 10.6, 21.8, 32.5, 43.1];
n_concs = length(xlink_concs);
norm_intensity = zeros([n_datapoints n_concs]);
norm_intensity_naught = zeros([n_datapoints n_concs]);

for i_run=1:1:2
for i_conc=1:1:n_concs
    if i_run == 1
        simName = sprintf('coop_HiKD_HiC_0_%i', int32(xlink_concs(i_conc)*10));    
    else 
        simName = sprintf('coop_HiKD_HiC_%i', int32(xlink_concs(i_conc)*10));
    end
    
    fileName = sprintf(fileDirectory, sprintf(fileStructure, simName));
    xlink_data_file = fopen(fileName);
    xlink_raw_data = fread(xlink_data_file, [n_mts * n_sites_max * n_datapoints], '*int');
    fclose(xlink_data_file);
    xlink_data = reshape(xlink_raw_data, n_sites_max, n_mts, n_datapoints);
    
    xlink_data(xlink_data ~= xlink_speciesID) = 0;
    xlink_data(xlink_data == xlink_speciesID) = 1;
    
    for i_data=1:1:n_datapoints
        for i_mt=1:1:n_mts
            occupancy = xlink_data(:, i_mt, i_data);
            n_sites = mt_length(i_mt);
            for i_site=1:1:n_sites
                if i_run == 1
                norm_intensity_naught(i_data, i_conc) = norm_intensity_naught(i_data, i_conc) + occupancy(i_site);
                else
                norm_intensity(i_data, i_conc) = norm_intensity(i_data, i_conc) + occupancy(i_site);    
                end
            end
        end
    end
end
end

norm_intensity = smoothdata(norm_intensity);
intensity_avgs = zeros([n_concs 1]);
intensity_avgs_naught = zeros([n_concs 1]);
for i_conc=1:1:n_concs
   intensity_avgs(i_conc) =  mean(norm_intensity(5000:10000, i_conc));
   intensity_avgs_naught(i_conc) = mean(norm_intensity_naught(5000:10000, i_conc));
end

exp_prediction = intensity_avgs(1) * exp_scaling;
plot(xlink_concs, exp_prediction, '^', 'LineWidth', 2, 'MarkerSize', 14);
hold on
plot(xlink_concs, intensity_avgs, '*', 'LineWidth', 2, 'MarkerSize', 14);
plot(xlink_concs, intensity_avgs_naught, '*', 'LineWidth', 2, 'MarkerSize', 14);

%title('Interaction energy of -2.5 kBT for PRC1-PRC1 neighbors');
xlabel('PRC1 Concentration (nM)', 'FontSize', 14);
ylabel('Intensity (A.U.)', 'FontSize', 14);
legend({'Experiment', 'Simulation (-2.75 kbT)', 'Simulation (0 kbT)'}, 'location', 'northwest', 'FontSize', 14);