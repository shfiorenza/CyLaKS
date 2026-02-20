clear variables;

sim_name = 'test2b';

% Load parameter structure
file_dir = '..';  % Default; only change if you move CyLaKS output files
params = load_parameters(sprintf('%s/%s', file_dir, sim_name));

filament_num_file = fopen(sprintf('%s/%s_filament_num.file', file_dir, sim_name));
filament_num_raw = fread(filament_num_file, '*double');
fclose(filament_num_file);
filament_forces_file = fopen(sprintf('%s/%s_filament_forces.file', file_dir, sim_name));
filament_forces_raw = fread(filament_forces_file, '*double');
fclose(filament_forces_file);

mt_num = filament_num_raw; 
mt_forces_x = NaN(params.n_datapoints, max(mt_num));
mt_forces_y = NaN(params.n_datapoints, max(mt_num));

i_data = 1;
for i_datapoint = 1 : 1 : params.n_datapoints
    for i_mt = 1 : 1 : mt_num(i_datapoint)
        mt_forces_x(i_datapoint, i_mt) = filament_forces_raw(i_data);
        i_data = i_data+1;
        mt_forces_y(i_datapoint, i_mt) = filament_forces_raw(i_data);
        i_data = i_data+1;
    end
end


% Open figure and set to desired size (each frame must be this same size)
fig1 = figure('Position', [50 50 1000 500]);
plot(linspace(0, params.t_run, params.n_datapoints), mt_forces_x(:, 1));
hold on
for i_mt = 2 : max(mt_num)
    plot(linspace(0, params.t_run, params.n_datapoints), mt_forces_x(:, i_mt))
end
plot(linspace(0, params.t_run, params.n_datapoints), mean(mt_forces_x, 2, "omitnan"), 'LineWidth', 3, 'Color', 'black')

fig2 = figure('Position', [150 50 1000 500]);
plot(linspace(0, params.t_run, params.n_datapoints), mt_forces_y(:, 1));
hold on
for i_mt = 2 : max(mt_num)
    plot(linspace(0, params.t_run, params.n_datapoints), mt_forces_y(:, i_mt))
end


%fig = figure();
%plot(mt_forces_y, params.n_datapoints);
