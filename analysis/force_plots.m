clear variables;
% *** modify these variables ***
sim_name = 'demo_75_5';
start_time = 0;  % in seconds 
end_time = -1;   % in seconds; -1 defaults to full data 
smoothing_window  = 1;   % in seconds

% *** automatically assigned variables -- do not modify ***
% Open log file and parse it into param labels & their values
fileDirectory = '../%s';
%fileDirectory = '/home/shane/data_kif4a_paper/run_mobility_both/%s';
log_file = sprintf(fileDirectory, sprintf('%s.log', sim_name));
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
n_dims = 2; % hard-coded for now; CyLaKS always outputs data in 2-D


xlinkFileName = '%s_xlink_force.file';
xlink_data_file = fopen(sprintf(fileDirectory, sprintf(xlinkFileName, sim_name)));
raw_xlink_data = fread(xlink_data_file, [n_dims, n_datapoints], 'double');
fclose(xlink_data_file);

if(end_time == -1)
    end_time = n_datapoints * time_per_datapoint;
end
i_start = int32(start_time / time_per_datapoint) + 1;
i_end = int32(end_time / time_per_datapoint);
n_datapoints_active = i_end - i_start + 1;

smoothing_window_size = smoothing_window / time_per_datapoint;

xlink_data = raw_xlink_data; % (:, i_start:i_end);
smoothed_x = smooth(xlink_data(1, :), smoothing_window_size);
smoothed_y = smooth(xlink_data(2, :), smoothing_window_size);
smoothed_tot = smooth(sqrt(xlink_data(1, :).*xlink_data(2, :)), smoothing_window_size);

fig1 = figure();
set(fig1, 'Position', [100, 100, 1000, 500]);
hold on
plot(linspace(start_time, end_time, n_datapoints_active), smoothed_x(i_start:i_end), 'LineWidth', 2);
plot(linspace(start_time, end_time, n_datapoints_active), smoothed_y(i_start:i_end), 'LineWidth', 2);
%plot(linspace(start_time, end_time, n_datapoints), smoothed_tot, 'LineWidth', 2);


ylabel('Applied force (pN)');
xlabel('Time (s)');
% cip off boundaries b/c smoothing makes them artificially large 
xlim([start_time+smoothing_window end_time-smoothing_window]);
grid off
legend({'Horizontal component', 'Vertical component'}, 'location', 'best', 'FontSize', 12);
legend('boxoff');