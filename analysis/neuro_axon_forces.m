clear variables;

sim_name = 'test';

% Load parameter structure
file_dir = '..';  % Default; only change if you move CyLaKS output files
params = load_parameters(sprintf('%s/%s', file_dir, sim_name));


filament_filename = sprintf('%s/%s_axon_coords.file', file_dir, sim_name);
filament_file = fopen(filament_filename);
data_raw = fread(filament_file, '*double');
n_mts = zeros(1, params.n_datapoints);

forces_filename = sprintf('%s/%s_axon_forces.file', file_dir, sim_name);
forces_file = fopen(forces_filename);
data_forces_raw = fread(forces_file, '*double');

i_data = 1;
for i_datapoint = 1 : 1 : params.n_datapoints
    n_mts(i_datapoint) = data_raw(i_data);
    i_data = i_data + 1;
    for i_mt = 1 : 1 : n_mts(i_datapoint)
        i_data = i_data+1;
        for i_dim = 1 : 1 : params.n_dims
            for i_end = 1 : 1 : 2
                i_data = i_data + 1;
            end
        end
    end
end

i_data = 1;
for i_datapoint = 1 : 1 : params.n_datapoints
    for i_mt = 1 : 1 : n_mts(i_datapoint)
        mt_forces_x(i_mt, i_datapoint) = data_forces_raw(i_data);
        i_data = i_data+1;
        mt_forces_y(i_mt, i_datapoint) = data_forces_raw(i_data);
        i_data = i_data+1;
    end
    for i_mt = n_mts(i_datapoint) + 1 : max(n_mts)
        mt_forces_x(i_mt, i_datapoint) = nan;
        mt_forces_y(i_mt, i_datapoint) = nan;
    end
end


% Open figure and set to desired size (each frame must be this same size)
fig1 = figure('Position', [50 50 1000 500]);
plot(linspace(0, params.t_run, params.n_datapoints), mt_forces_x(1, :));
hold on
for i_mt = 2 : max(n_mts)
    plot(linspace(0, params.t_run, params.n_datapoints), mt_forces_x(i_mt, :))
end
plot(linspace(0, params.t_run, params.n_datapoints), mean(mt_forces_x, 1, "omitnan"), 'LineWidth', 3, 'Color', 'black')

fig2 = figure('Position', [150 50 1000 500]);
plot(linspace(0, params.t_run, params.n_datapoints), mt_forces_y(1, :));
hold on
for i_mt = 2 : max(n_mts)
    plot(linspace(0, params.t_run, params.n_datapoints), mt_forces_y(i_mt, :))
end


%fig = figure();
%plot(mt_forces_y, params.n_datapoints);
