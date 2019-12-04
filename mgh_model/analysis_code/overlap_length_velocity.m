clear all;
% Input parameters
simName = "2019_12_02_slideScan/slide_scan_600_6";
start_point = 0;
site_size = 0.008; % in microns
% Set file directory & generate .log file name
fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
log_file = sprintf(fileDirectory, sprintf('%s.log', simName));
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
% Read in system parameters
n_steps = str2double(values{contains(params, "n_steps")});
delta_t = sscanf(values{contains(params, "delta_t")}, '%g');
n_datapoints = str2double(values{contains(params, "n_datapoints")});
% Calculate parameters for plotting / etc;
time_per_datapoint = delta_t * n_steps / n_datapoints;
end_time = n_datapoints * time_per_datapoint;
start_time = start_point * time_per_datapoint;
% Open mt coordinate file
mt_coords_file = fopen(sprintf(fileDirectory, sprintf('%s_mt_coord.file', simName)));
mt_coord_data = fread(mt_coords_file, [n_mts, n_datapoints], 'double');
fclose(mt_coords_file);
% Run through raw coord data to get overlap length at every datapoint
final_overlap_data = zeros(n_datapoints, 1);
end_dist_data = zeros(n_datapoints, 1);
for i_data = start_point + 1: n_datapoints
    mt_coords = mt_coord_data(:, i_data)';
    mt_endpoints = mt_coords + mt_lengths;
    %{
    overlap_start = max(mt_coords);
    overlap_end = min(mt_endpoints);
    overlap_length = (overlap_end - overlap_start) * site_size;
    %}
    plus_end_one = mt_coords(1);
    plus_end_two = mt_coords(2) + mt_lengths(2);
  
    final_overlap_data(i_data) = abs(plus_end_two - plus_end_one) * site_size; 
    
end

% Use gradient function with above spacing to get slope of overlap length
final_overlap_data = smoothdata(final_overlap_data, 'movmean', 10);
slope_data = smoothdata(gradient(final_overlap_data, time_per_datapoint), 'movmean', 100);

%fig1 = figure();
%set(fig1, 'Position', [50, 50, 2.5*480, 2.5*300])
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
ylim([-100 100]);
%axis tight
grid on
grid minor

set(gcf, 'color', 'w');