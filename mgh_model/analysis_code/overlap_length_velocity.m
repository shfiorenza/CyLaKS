clear all;
% Often-changed variables
mt_lengths = [1000, 500];
simName = '2019_11_14_slideScan/slide_scan_500_2500_5';
%simName = 'test_bias'
% Pseudo-constant variables
n_steps = 60000000;
delta_t = 0.000025;
n_datapoints = 10000;
start_point = 0;
site_size = 0.008; % in microns
% Calculate parameters for plotting / etc;
n_mts = length(mt_lengths);
time_per_datapoint = delta_t * n_steps / n_datapoints;
end_time = n_datapoints * time_per_datapoint;
start_time = start_point * time_per_datapoint;

fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
mtCoordFileName = '%s_mt_coord.file';

mt_coords_file = fopen(sprintf(fileDirectory, sprintf(mtCoordFileName, simName)));
mt_coord_data = fread(mt_coords_file, [n_mts, n_datapoints], 'double');
%mt_coord_data = reshape(mt_coord_data, n_datapoints, n_mts);
fclose(mt_coords_file);

% Run through raw coord data to get overlap length at every datapoint
final_overlap_data = zeros(n_datapoints, 1);
end_dist_data = zeros(n_datapoints, 1);
for i_data = start_point + 1: n_datapoints
    mt_coords = mt_coord_data(:, i_data)';
    mt_endpoints = mt_coords + mt_lengths;
    
    overlap_start = max(mt_coords);
    overlap_end = min(mt_endpoints);
    
    overlap_length = (overlap_end - overlap_start) * site_size;
  
    final_overlap_data(i_data) = overlap_length; 
    
    %{
    plus_end_one = mt_coord_one;
    plus_end_two = mt_coord_two + n_sites(2);
    end_dist = (plus_end_two - plus_end_one) * 0.008;
    end_dist_data(i_data) = end_dist;
    %}
end

% Use gradient function with above spacing to get slope of overlap length
final_overlap_data = smoothdata(final_overlap_data, 'movmean', 10);
slope_data = smoothdata(gradient(final_overlap_data, time_per_datapoint), 'movmean', 100);

fig1 = figure();
set(fig1, 'Position', [50, 50, 2.5*480, 2.5*300])
hold on

% Plot overlap length on top
subplot(2, 1, 1)
hold all
plot(linspace(start_time, end_time, n_datapoints), final_overlap_data, 'LineWidth', 2);
%title(sprintf('Overlap length over time (%g microns or %d sites in length)', ...
 %length, n_sites(2)));
 %title('c eff teth = 1000 nM; k on xlink = 0.1 (nM * s)^{-1}');
ylabel('Overlap length (microns)');
xlabel('Time (s)');
axis tight
xlim([start_time end_time]);
ylim([0 max(mt_lengths) * 0.008]);
grid on
grid minor

% Plot sliding velocity on bottom
subplot(2, 1, 2)
plot(linspace(start_time, end_time, n_datapoints), slope_data * 1000, ...
        'LineWidth', 2);
title('Sliding velocity over time');
ylabel('Sliding velocity (nm/s)');
xlabel('Time (s)');
xlim([start_time end_time]);
%ylim([-250 250]);
%axis tight
grid on
grid minor

set(gcf, 'color', 'w');