clear variables;

sim_name = 'test_longID';

output_movie_name = 'test_long2_profile';

start_frame = 1;
end_frame = -1;  % set to -1 to run until end of data

frames_per_plot = 10;
movie_duration = 15; % in seconds

% Load parameter structure
file_dir = '..';  % Default; only change if you move CyLaKS output files
params = load_parameters(sprintf('%s/%s', file_dir, sim_name));

if end_frame == -1
    end_frame = params.n_datapoints;
end
active_frames = end_frame - start_frame;

% Initialize videowriter object
v = VideoWriter(output_movie_name, 'MPEG-4');
v.FrameRate = (active_frames / frames_per_plot) / movie_duration;
open(v);
frame_box = [0 0 1445 200];


filament_num_file = fopen(sprintf('%s/%s_filament_num.file', file_dir, sim_name));
filament_num_raw = fread(filament_num_file, '*double');
fclose(filament_num_file);
filament_lengths_file = fopen(sprintf('%s/%s_filament_lengths.file', file_dir, sim_name));
filament_lengths_raw = fread(filament_lengths_file, '*double');
fclose(filament_lengths_file);
filament_pos_file = fopen(sprintf('%s/%s_filament_pos.file', file_dir, sim_name));
filament_pos_raw = fread(filament_pos_file, '*double');
fclose(filament_pos_file);
filament_tip_file = fopen(sprintf('%s/%s_filament_pos_tip.file', file_dir, sim_name));
filament_tip_raw = fread(filament_tip_file, '*double');
fclose(filament_tip_file);

mt_num = filament_num_raw;
mt_len = NaN(params.n_datapoints, max(mt_num));
mt_pos = NaN(params.n_datapoints, max(mt_num), 2, params.n_dims);
mt_tip = NaN(params.n_datapoints, max(mt_num), params.n_dims);


i_entry = 1;
j_entry = 1;
k_entry = 1;
for i_datapoint = 1 : 1 : params.n_datapoints
    for i_mt = 1 : mt_num(i_datapoint)
        mt_len(i_datapoint, i_mt) = filament_lengths_raw(i_entry);
        i_entry = i_entry + 1;
        for i_end = 1 : 1 : 2
            for i_dim = 1 : 1 : params.n_dims
                mt_pos(i_datapoint, i_mt, i_end, i_dim) = filament_pos_raw(j_entry);
                j_entry = j_entry + 1;
            end
        end
        for i_dim = 1 : 1 : params.n_dims
            mt_tip(i_datapoint, i_mt, i_dim) = filament_tip_raw(k_entry);
            k_entry = k_entry + 1;
        end
    end
end

% Open figure and set to desired size (each frame must be this same size)
fig1 = figure('Position', [50 50 1000 250]);


min_x = min(mt_pos(:, :, :, 1), [], "all");
max_x = max(mt_pos(:, :, :, 1), [], "all");
span = max_x - min_x;
n_sites = span / 8.2;
occu_data_stable = zeros(active_frames, ceil(n_sites));
occu_data_unstable = zeros(active_frames, ceil(n_sites));
% Run through all datapoints; each one is a frame in our movie
for i_frame = start_frame : frames_per_plot : end_frame - 1
    % Clear figure so that it only displays figures from current datapoint
    clf;
    % Set Axes properties
    ax = axes('Units', 'normalized', 'Position', [0.075 0.17 0.9 0.8]);
    set(gca,'xdir','reverse');%,'ydir','reverse')
    hold all;
    for i_data = i_frame : i_frame + frames_per_plot - 1
    %i_data = i_frame;
    for i_mt = 1 : mt_num(i_data)
        plus_pos = mt_pos(i_data, i_mt, 1, :);
        minus_pos = mt_pos(i_data, i_mt, 2, :);
        tip_pos = mt_tip(i_data, i_mt, :);
        len_sites = mt_len(i_data, i_mt);
        %disp(len_sites * 8.2)
        %l_sites = abs(plus_pos(1) - minus_pos(1))
        len_tip = abs(tip_pos(1) - plus_pos(1));
        pos_start = min(mt_pos(i_data, i_mt, :, 1));
        i_start = ceil((pos_start - min_x)/8.2) + 1;
        for i = i_start : i_start + len_sites - 2
            if len_tip == 0
                occu_data_stable(i_frame, i) = occu_data_stable(i_frame, i) + 1/frames_per_plot;
            else
                occu_data_unstable(i_frame, i) = occu_data_unstable(i_frame, i) + 1/frames_per_plot;
            end
        end
    end
    end
    plot(linspace(min_x/1000.0, max_x/1000.0, ceil(n_sites)), occu_data_stable(i_frame, :), 'LineWidth', 3)
    hold on
    plot(linspace(min_x/1000.0, max_x/1000.0, ceil(n_sites)), occu_data_unstable(i_frame, :), 'LineWidth', 3)
    legendLabel = {"Stable MTs", "Unstable MTs"};
    legend(legendLabel, 'location', 'northeast');
    ylim([0 300]);
    xlim([-6 16]);
    xlabel("Position (microns)");
    ylabel("Microtubule density (A.U.)")
    xlim([min_x/1000.0 - 1 max_x/1000.0 + 1])
    dim = [0.11 0.625 .3 .3];
    time = (i_frame - start_frame) * params.time_per_datapoint;
    str = sprintf('Time: %#.2f seconds', time);
    annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on');
    drawnow();
    writeVideo(v, getframe(gcf));
end

close(v);
