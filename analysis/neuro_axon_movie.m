clear variables;

sim_name = 'test_long';

output_movie_name = 'test_long2';

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
r_prot = (params.site_size*1000);

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

% Run through all datapoints; each one is a frame in our movie
for i_data = start_frame : frames_per_plot : end_frame
    % Clear figure so that it only displays figures from current datapoint
    clf;
    % Set Axes properties
    ax = axes('Units', 'normalized', 'Position', [0.075 0.15 0.9 0.8]);
    set(gca,'xdir','reverse');%,'ydir','reverse')
    hold all;
    min_x = min(min(mt_pos(i_data, :, :, 1)));
    max_x = max(max(mt_pos(i_data, :, :, 1)));
    min_y = min(min(mt_pos(i_data, :, :, 2)));
    max_y = max(max(mt_pos(i_data, :, :, 2)));
    %ax.XLim = [(min_x - 250) (max_x + 250)];
    %ax.XLim = [-5500 10500];
    %ax.XLim = [-7500 45000];
    %ax.XLim = [-100000 100000];
    %ax.YLim = [(min_y - 500) (max_y + 500)];
    ax.YLim = [-40 100];
    ax.XLim = [-1000 41000];
    %ax.TickLength = [0.005 0.005];
    ax.XLabel.String = 'x position (nm)';
    ax.YLabel.String = 'y position (nm)';
    % Draw filaments
    % if(params.n_mts > 1)
    %     com_y_one = (filament_pos(2, 1, 1, i_data) + filament_pos(2, 2, 1, i_data))/2;
    %     com_y_two = (filament_pos(2, 1, 2, i_data) + filament_pos(2, 2, 2, i_data))/2;
    % end
    for i_mt = 1:1:mt_num(i_data)
        plus_pos = mt_pos(i_data, i_mt, 1, :);
        minus_pos = mt_pos(i_data, i_mt, 2, :);
        tip_pos = mt_tip(i_data, i_mt, :);
        if plus_pos(1) > minus_pos(1)
            polarity = 0;
            color = [0.7 0.7 0.7];
        else
            polarity = 1;
            color = [0.25 0.25 0.25];
        end
        %plus_pos = filament_pos(:, 1, i_mt, i_data);
        %minus_pos = filament_pos(:, 2, i_mt, i_data);
        line([plus_pos(1)-r_prot/2, minus_pos(1)-r_prot/2],[plus_pos(2), minus_pos(2)], ...
            'LineWidth', 2, 'Color', color);
        line([plus_pos(1)-r_prot/2, tip_pos(1)-r_prot/2],[plus_pos(2), tip_pos(2)], ...
            'LineWidth', 2, 'Color', [0 1 0]);
        %rectangle('Position', [plus_pos(1)-r_prot/2 plus_pos(2)-2*r_prot r_prot 4*r_prot], ...
        %     'FaceColor', [0 0 0], 'Curvature', [1 1]);
        %n_sites = params.mt_lengths(i_mt);
        %dx = -1;
        %mt_dir = 1;
        % line_vec = [minus_pos(1) - plus_pos(1), minus_pos(2) - plus_pos(2)];
        % if params.polarity(i_mt) == 1
        %     dx = 1;
        %     mt_dir = -1;
        %     line_vec = [plus_pos(1) - minus_pos(1), plus_pos(2) - minus_pos(2)];
        % end
    end
    dim = [0.085 0.645 .3 .3];
    time = (i_data - start_frame) * params.time_per_datapoint;
    str = sprintf('Time: %#.2f seconds', time);
    annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on');
    drawnow();
   writeVideo(v, getframe(gcf));
end

close(v);
