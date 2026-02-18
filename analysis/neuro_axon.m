clear variables;

sim_name = 'test';

output_movie_name = 'test';

start_frame = 1; 
end_frame = -1;  % set to -1 to run until end of data

frames_per_plot = 100; 
movie_duration = 30; % in seconds

% Load parameter structure
file_dir = '..';  % Default; only change if you move CyLaKS output files
params = load_parameters(sprintf('%s/%s', file_dir, sim_name));

if end_frame == -1
    end_frame = params.n_datapoints; 
end
active_frames = end_frame - start_frame;
r_prot = (params.site_size*1000);

% Initialize videowriter object
v = VideoWriter(output_movie_name);%, 'MPEG-4');
v.FrameRate = (active_frames / frames_per_plot) / movie_duration;
open(v);
frame_box = [0 0 1445 200];


% Open figure and set to desired size (each frame must be this same size)
fig1 = figure('Position', [50 50 1000 500]);

filament_filename = sprintf('%s/%s_axon_coords.file', file_dir, sim_name);
filament_file = fopen(filament_filename);
data_raw = fread(filament_file, '*double');
n_mts = zeros(1, params.n_datapoints);
i_data = 1;
for i_datapoint = 1 : 1 : params.n_datapoints
    n_mts(i_datapoint) = data_raw(i_data);
    i_data = i_data + 1;
    for i_mt = 1 : 1 : n_mts(i_datapoint)
        mt_lengths(i_mt, i_datapoint) = data_raw(i_data);
        i_data = i_data+1;
        for i_dim = 1 : 1 : params.n_dims
            for i_end = 1 : 1 : 2
                filament_pos(i_end, i_dim, i_mt, i_datapoint) = data_raw(i_data);
                i_data = i_data + 1;
            end
        end
    end
end


% Run through all datapoints; each one is a frame in our movie
for i_data = start_frame : frames_per_plot : end_frame
    % Clear figure so that it only displays figures from current datapoint
    clf;
    % Set Axes properties
    ax = axes('Units', 'normalized', 'Position', [0.075 0.085 0.9 0.9]);
    set(gca,'xdir','reverse');%,'ydir','reverse')
    hold all;
    min_x = min(min(filament_pos(1, :, :, i_data)));
    max_x = max(max(filament_pos(1, :, :, i_data)));
    min_y = min(min(filament_pos(2, :, :, i_data)));
    max_y = max(max(filament_pos(2, :, :, i_data)));
    %ax.XLim = [(min_x - 25) (max_x + 25)];
    ax.XLim = [-7500 200000];
    %ax.XLim = [-100000 100000];
    %ax.YLim = [(min_y - 500) (max_y + 500)];
    ax.YLim = [-1000 1000];
    ax.TickLength = [0.005 0.005];
    ax.XLabel.String = 'x position (nm)';
    ax.YLabel.String = 'y position (nm)';
    % Draw filaments
    % if(params.n_mts > 1)
    %     com_y_one = (filament_pos(2, 1, 1, i_data) + filament_pos(2, 2, 1, i_data))/2;
    %     com_y_two = (filament_pos(2, 1, 2, i_data) + filament_pos(2, 2, 2, i_data))/2;
    % end
    for i_mt = 1:1:n_mts(i_data)
        plus_pos = filament_pos(:, 1, i_mt, i_data);
        minus_pos = filament_pos(:, 2, i_mt, i_data);
        line([plus_pos(1)-r_prot/2, minus_pos(1)-r_prot/2],[plus_pos(2), minus_pos(2)], ...
            'LineWidth', 2, 'Color', [0.7 0.7 0.7]);
        rectangle('Position', [plus_pos(1)-r_prot/2 plus_pos(2)-2*r_prot r_prot 4*r_prot], ...
             'FaceColor', [0 0 0], 'Curvature', [1 1]);
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
    dim = [0.11 0.625 .3 .3];
    time = (i_data - start_frame) * params.time_per_datapoint;
    str = sprintf('Time: %#.2f seconds', time);
    annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on');
    drawnow();
   writeVideo(v, getframe(gcf));
end

close(v);
