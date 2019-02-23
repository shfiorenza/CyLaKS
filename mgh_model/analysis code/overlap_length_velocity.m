clear all;
% Often-changed variables
n_sites = 500;
simName = 'slide_scan1/unteth_0.005';
% Pseudo-constant variables
n_mts = 2;
n_steps = 100000000;
delta_t = 0.00001;
n_datapoints = 10000;
starting_point = 0;
% Calculate parameters for plotting / etc;
length = n_sites * 0.008;
end_time = n_steps * delta_t;
start_time = starting_point * delta_t;

%fileDirectory = '/home/shane/Desktop/slide_scan/%s';
fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
fileStructure = '%s_mt_coord.file';
fileName = sprintf(fileDirectory, sprintf(fileStructure, simName));

data_file = fopen(fileName);
raw_data = fread(data_file, [n_mts, n_datapoints], 'double');
fclose(data_file);

% Run through raw coord data to get overlap length at every datapoint
final_overlap_data = zeros(n_datapoints, 1);
for i_data = starting_point + 1 : 1 : n_datapoints
    mt_coord_one = raw_data(1, i_data);
    mt_coord_two = raw_data(2, i_data);
    delta = (mt_coord_two - mt_coord_one) * 0.008;
    if(delta > 0)
        overlap_length = length - delta;
    else
        overlap_length = length + delta;
    end
    final_overlap_data(i_data) = overlap_length; 
end

% Calculate real time that passes per iteration to get an accurate velocity
time_per_datapoint = delta_t * (n_steps / n_datapoints);
% Use gradient function with above spacing to get slope of overlap length
smoothed_overlap_data = smooth(final_overlap_data, 10);
%smoothed_overlap_data = final_overlap_data;
slope_data = gradient(smoothed_overlap_data, time_per_datapoint);
final_slope_data = smooth(abs(slope_data), 100);
final_overlap_data = smoothed_overlap_data;

fig1 = figure();
set(fig1, 'Position', [50, 50, 2.5*480, 2.5*300])

% Plot overlap length on top
subplot(2, 1, 1)
plot(linspace(start_time, end_time, n_datapoints), final_overlap_data, ...
        'LineWidth', 2);
title(sprintf('Overlap length over time (%g microns or %d sites in length)', ...
        length, n_sites));
ylabel('Overlap length (microns)');
xlabel('Time (s)');
axis tight
xlim([start_time + 1 end_time - 1]);
ylim([0 n_sites * 0.008]);
grid on
grid minor

% Plot sliding velocity on bottom
subplot(2, 1, 2)
plot(linspace(start_time, end_time, n_datapoints), final_slope_data * 1000, ...
        'LineWidth', 2);
title('Sliding velocity over time');
ylabel('Sliding velocity (nm/s)');
xlabel('Time (s)');
xlim([start_time + 1 end_time - 1]);
%ylim([-500 500]);
%axis tight
grid on
grid minor

set(gcf, 'color', 'w');