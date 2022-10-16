clear all;
% Input parameters
sim_name = "shep_slide_1nM_200nM_duo_1250";
seeds = [0 1 2 3]; % 3]; %,2,5];%,6,7,4,8];%,8,9]%,11];
start_point = 0;


% Set file directory
%fileDirectory = '/home/shane/Projects/CyLaKS/%s';
fileDirectory = '../%s';
log_file = sprintf(fileDirectory, sprintf('%s.log', sim_name + "_0"));
log = textscan(fileread(log_file), '%s %s', 'Delimiter', '=');
params = log{1, 1};
values = log{1, 2};
% Read in system params
dt = sscanf(values{contains(params, "dt ")}, '%g');
time_per_datapoint = sscanf(values{contains(params, "t_snapshot ")}, '%g');
n_datapoints = str2double(values{contains(params, "n_datapoints ")});
% Use actual recorded number of datapoints to parse thru data/etc
if any(contains(params, "N_DATAPOINTS ") ~= 0)
    n_datapoints = str2double(values{contains(params, "N_DATAPOINTS ")});
end
site_size =  sscanf(values{contains(params, "site_size ")}, '%g') / 1000; % in um
% Read in number of MTs
n_mts = sscanf(values{contains(params, "count ")}, '%g');
if any(contains(params, "COUNT ") ~= 0)
    n_mts = sscanf(values{contains(params, "COUNT ")}, '%g');
end
% Read in MT lengths (in n_sites)
mt_lengths = zeros(1, n_mts);
for i_mt = 1 : n_mts
    string = sprintf("n_sites[%i] ", i_mt - 1);
    mt_lengths(i_mt) = sscanf(values{contains(params, string)}, '%i');
    if any(contains(params, sprintf("N_SITES[%i] ", i_mt - 1)) ~= 0)
        string = sprintf("N_SITES[%i] ", i_mt - 1);
        mt_lengths(i_mt) = sscanf(values{contains(params, string)}, '%i');
    end
end
max_sites = max(mt_lengths);
n_dims = 2; % hard-coded for now; CyLaKS always outputs data in 2-D
end_time = n_datapoints * time_per_datapoint;
start_time = start_point * time_per_datapoint;
active_datapoints = n_datapoints - start_point;

overlap_length_data = zeros(length(seeds), n_datapoints);
overlap_velocity = zeros(length(seeds), n_datapoints);
plus_end_dist_data = zeros(length(seeds), n_datapoints);
plus_end_velocity = zeros(length(seeds), n_datapoints);

% Run through raw coord data to get overlap length at every datapoint
for i_seed = 1:length(seeds)
    simName = sim_name;

    if length(seeds) > 1 % && i_seed > 1
        simName = sprintf("%s_%i", sim_name, seeds(i_seed));
    end

    % Open mt coordinate file
    filamentFileName = '%s_filament_pos.file';
    filamentFile = sprintf(fileDirectory, sprintf(filamentFileName, simName));
    filament_pos = zeros(n_dims, 2, n_mts, n_datapoints);
    if isfile(filamentFile)
        file = fopen(filamentFile);
        data = fread(file, 2*n_dims * n_mts * n_datapoints, '*double');
        fclose(file);
        filament_pos = reshape(data, n_dims, 2, n_mts, n_datapoints);
    end
    
    for i_data = start_point + 1:n_datapoints
        
       % plus_pos = filament_pos(:, 1, i_mt, i_data);
        %minus_pos = filament_pos(:, 2, i_mt, i_data);
        
        %mt_coords = filament_pos(:, i_data)';
        %mt_endpoints = mt_coords + mt_lengths;
        %overlap_start = max(mt_coords);
        %overlap_end = min(mt_endpoints);
        %overlap_length = (overlap_end - overlap_start) * site_size;
        %overlap_length_data(i_seed, i_data) = overlap_length;

        plus_end_one = filament_pos(1, 1, 1, i_data);
        plus_end_two = filament_pos(1, 1, 2, i_data);
        plus_end_dist = abs(plus_end_two - plus_end_one);
        plus_end_dist_data(i_seed, i_data) = plus_end_dist;

    end

end

% Use gradient function with above spacing to get slope of overlap length
for i_seed = 1:length(seeds)
    %overlap_length_data(i_seed, :) = smoothdata(overlap_length_data(i_seed, :), 'movmean', 10);
    %overlap_velocity(i_seed, :) = smoothdata(gradient(overlap_length_data(i_seed, :), ...
    %    time_per_datapoint), 'movmean', 100);
    % Convert velocity from um/s to nm/s
    %overlap_velocity(i_seed, :) = overlap_velocity(i_seed, :) * 1000;

    plus_end_dist_data(i_seed, :) = smoothdata(plus_end_dist_data(i_seed, :), 'movmean', 500);
    plus_end_velocity(i_seed, :) = smoothdata(gradient(plus_end_dist_data(i_seed, :), ...
        time_per_datapoint), 'movmean', 500);
    % Convert velocity from um/s to nm/s
    plus_end_velocity(i_seed, :) = plus_end_velocity(i_seed, :);
end
fig1 = figure();
set(fig1, 'Position', [50, 50, 2.5 * 300, 2 * 300])
hold all


% Plot overlap length on top
subplot(2, 1, 1)
%plot(linspace(start_time, end_time, active_datapoints), overlap_length_data, 'LineWidth', 2);
hold on
plot(linspace(start_time, end_time, active_datapoints), plus_end_dist_data * 0.001, 'LineWidth', 2);

%title('Distance between plus-ends over time');
ylabel('Plus-end distance (\mum)');
%xlabel('Time (s)');
axis tight
xlim([start_time end_time]);
ylim([0 (11/10)*max(max(plus_end_dist_data* 0.001))]);
set(gca, 'FontSize', 14);


% Plot sliding velocity on bottom
subplot(2, 1, 2)
%plot(linspace(start_time, end_time, active_datapoints), overlap_velocity, 'LineWidth', 2);
hold on
plot(linspace(start_time, end_time, active_datapoints), plus_end_velocity * 0.001, 'LineWidth', 2);
%plot(linspace(start_time, end_time, active_datapoints), gradient(plus_end_velocity * 0.001, 0.008), 'LineWidth', 2);
%title('Sliding velocity over time');
ylabel('Sliding velocity (\mum/s)');
xlabel('Time (s)');
xlim([start_time end_time]);
ylim([(11/10)*min(min(plus_end_velocity* 0.001)) (11/10)*max(max(plus_end_velocity* 0.001))]);
%axis tight

disp(mean(plus_end_dist_data(:, n_datapoints)));
%disp(mean(plus_end_velocity(10500:12000)));

%sgtitle('10 nM PRC1 + 2000 nM K401');
set(gca, 'FontSize', 14);

set(gcf, 'color', 'w');
