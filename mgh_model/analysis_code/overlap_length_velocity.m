clear all;
% Input parameters
baseName = "test_expandA";
seeds = [0]; %,2,5];%,6,7,4,8];%,8,9]%,11];
start_point = 0;
site_size = 0.008; % in nm
% Set file directory
fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';

log_file = sprintf(fileDirectory, sprintf('%s.log', baseName)); %, seeds(1)));
% Open log file and parse it into param labels & their values
log = textscan(fileread(log_file),'%s %s', 'Delimiter', '=');
params = log{1,1};
values = log{1,2};
% Read in number of MTs
n_mts = str2double(values{contains(params, "count")});
% Read in length of each MTs

[length_one, length_two] = values{contains(params, "length")};
mt_lengths(1) = sscanf(length_one, '%i');
mt_lengths(2) = sscanf(length_two, '%i');

%{
[l_1, l_2, l_3, l_4, l_5] = values{contains(params, "length")};
mt_lengths(1) = sscanf(l_1, '%i');
mt_lengths(2) = sscanf(l_2, '%i');
mt_lengths(3) = sscanf(l_3, '%i');
mt_lengths(4) = sscanf(l_4, '%i');
mt_lengths(5) = sscanf(l_5, '%i');
%}
% Read in system parameters
n_steps = str2double(values{contains(params, "n_steps")});
delta_t = sscanf(values{contains(params, "delta_t")}, '%g');
n_datapoints = str2double(values{contains(params, "n_datapoints")});
% Calculate parameters for plotting / etc;
time_per_datapoint = delta_t * n_steps / n_datapoints;
end_time = n_datapoints * time_per_datapoint;
start_time = start_point * time_per_datapoint;
active_datapoints = n_datapoints - start_point;

overlap_length_data = zeros(length(seeds), n_datapoints);
overlap_velocity = zeros(length(seeds), n_datapoints);
plus_end_dist_data = zeros(length(seeds), n_datapoints);
plus_end_velocity = zeros(length(seeds), n_datapoints);

% Run through raw coord data to get overlap length at every datapoint
for i_seed = 1:length(seeds)
    simName = baseName;
    if length(seeds) > 1
        simName = sprintf("%s_%i", baseName, seeds(i_seed));
    end
    % Open mt coordinate file
    mt_coords_file = fopen(sprintf(fileDirectory, sprintf('%s_mt_coord.file', simName)));
    mt_coord_data = fread(mt_coords_file, [n_mts, n_datapoints], 'double');
    fclose(mt_coords_file);
    for i_data = start_point + 1: n_datapoints
        mt_coords = mt_coord_data(:, i_data)';
        mt_endpoints = mt_coords + mt_lengths;
        
        overlap_start = max(mt_coords);
        overlap_end = min(mt_endpoints);
        overlap_length = (overlap_end - overlap_start) * site_size;
        overlap_length_data(i_seed, i_data) = overlap_length;
        
        plus_end_one = mt_coords(1);
        plus_end_two = mt_coords(2) + mt_lengths(2);
        plus_end_dist = abs(plus_end_two - plus_end_one) * site_size;
        plus_end_dist_data(i_seed, i_data) = plus_end_dist;
        
    end
end

% Use gradient function with above spacing to get slope of overlap length
for i_seed=1:length(seeds)
    overlap_length_data(i_seed, :) = smoothdata(overlap_length_data(i_seed, :), 'movmean', 10);
    overlap_velocity(i_seed, :) = smoothdata(gradient(overlap_length_data(i_seed, :), ...
        time_per_datapoint), 'movmean', 100);
    % Convert velocity from um/s to nm/s
    overlap_velocity(i_seed, :) = overlap_velocity(i_seed, :) * 1000;
    
    plus_end_dist_data(i_seed, :) = smoothdata(plus_end_dist_data(i_seed, :), 'movmean', 10);
    plus_end_velocity(i_seed, :) = smoothdata(gradient(plus_end_dist_data(i_seed, :), ...
        time_per_datapoint), 'movmean', 100);
    % Convert velocity from um/s to nm/s
    plus_end_velocity(i_seed, :) = plus_end_velocity(i_seed, :) * 1000;
end

fig1 = figure();
set(fig1, 'Position', [50, 50, 2.5*480, 2.5*300])
hold all

% Plot overlap length on top
subplot(2, 1, 1)
plot(linspace(start_time, end_time, active_datapoints), overlap_length_data, 'LineWidth', 2);
hold on
plot(linspace(start_time, end_time, active_datapoints), plus_end_dist_data, 'LineWidth', 2);
title('Distance between plus-ends over time');
ylabel('Plus-end distance (microns)');
xlabel('Time (s)');
axis tight
xlim([start_time end_time]);
%ylim([0 4.5]);
legend(["Overlap length", "Plus-end distance"]);
grid on
grid minor

% Plot sliding velocity on bottom
subplot(2, 1, 2)
plot(linspace(start_time, end_time, active_datapoints), overlap_velocity, 'LineWidth', 2);
hold on
plot(linspace(start_time, end_time, active_datapoints), plus_end_velocity, 'LineWidth', 2);
title('Sliding velocity over time');
ylabel('Sliding velocity (nm/s)');
xlabel('Time (s)');
xlim([start_time end_time]);
%ylim([-60 40]);
%axis tight
grid on
grid minor

sgtitle('0.2 nM PRC1 + 6 nM Kif4A');

set(gcf, 'color', 'w');