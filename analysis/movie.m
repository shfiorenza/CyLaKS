clear variables;

fileDirectory = '/home/shane/projects/CyLaKS/%s';

sim_name = 'test';

movie_name = 'test';
start_frame = 1; 
frames_per_plot = 10;
movie_duration = 30; % in seconds

% Open log file and parse it into param labels & their values
log_file = sprintf(fileDirectory, sprintf('%s.log', sim_name));
log = textscan(fileread(log_file), '%s %s', 'Delimiter', '=');
params = log{1, 1};
values = log{1, 2};
% Read in number of MTs
n_mts = str2double(values{contains(params, "count ")});
mt_lengths = zeros(1, n_mts);
for i_mt = 1 : n_mts
   string = sprintf("n_sites[%i] ", i_mt - 1);
   mt_lengths(i_mt) = sscanf(values{contains(params, string)}, '%i');
end
% Read in system params
dt = sscanf(values{contains(params, "dt ")}, '%g');
steps_per_datapoint = str2double(values{contains(params, "n_steps_per_snapshot ")});
time_per_datapoint = dt * steps_per_datapoint;
n_datapoints = str2double(values{contains(params, "n_datapoints ")});
% Use actual recorded number of datapoints to parse thru data/etc
if any(contains(params, "N_DATAPOINTS ") ~= 0)
    n_datapoints = str2double(values{contains(params, "N_DATAPOINTS ")});
end
n_dims = 2;
site_size = 0.008; % in um
max_sites = max(mt_lengths);
xlink_cutoff = 5;
teth_cutoff = 19;

% Colors
blue = [30 144 255] / 255;
purple = [128 0 128] / 255;

% Figure parameters (i.e., how they appear)
site_height = 1;
site_width = 1;
end_frame = n_datapoints;
active_frames = end_frame - start_frame;

% Videowriter details
v = VideoWriter(movie_name);
v.FrameRate = (active_frames / frames_per_plot) / movie_duration;
open(v);
frame_box = [0 0 1445 200];

% Figure details
fig1 = figure;
set(fig1, 'Position', [0 100 400 400]);

% File info
filamentFileName = '%s_filament_pos.file';
proteinFileName = '%s_protein_id.file';
tethFileName = '%s_tether_anchor_pos.file';
motorHeadFileName = '%s_motor_head_trailing.file';
filamentFile = sprintf(fileDirectory, sprintf(filamentFileName, sim_name));
proteinFile = sprintf(fileDirectory, sprintf(proteinFileName, sim_name));
tethFile = sprintf(fileDirectory, sprintf(tethFileName, sim_name));
motorHeadFile = sprintf(fileDirectory, sprintf(motorHeadFileName, sim_name));

filament_pos = zeros(n_dims, 2, n_mts, n_datapoints);
if isfile(filamentFile)
    file = fopen(filamentFile);
    data = fread(file, 2*n_dims * n_mts * n_datapoints, '*double');
    fclose(file);
    filament_pos = reshape(data, n_dims, 2, n_mts, n_datapoints);
end
protein_ids = zeros(max_sites, n_mts, n_datapoints) - 1;
if isfile(proteinFile)
    file = fopen(proteinFile);
    data = fread(file, n_mts * max_sites * n_datapoints, '*int');
    fclose(file);
    protein_ids = reshape(data, max_sites, n_mts, n_datapoints);
end
teth_data = zeros(max_sites, n_mts, n_datapoints) - 1;
if isfile(tethFile)
    %{
    file = fopen(tethFile);
    data = fread(file, [n_mts * max_sites * n_datapoints], '*double');
    fclose(file);
    teth_data = reshape(data, max_sites, n_mts, n_datapoints);
    %}
end
motor_trailing = zeros(max_sites, n_mts, n_datapoints);
if isfile(motorHeadFile)
    %{
    file = fopen(motorHeadFile);
    data = fread(file, [n_mts * max_sites * n_datapoints], '*bool');
    fclose(file);
    motor_trailing = reshape(data, max_sites, n_mts, n_datapoints);
    %}
end

% Run through all datapoints; each one is a frame in our movie
for i_data = start_frame : frames_per_plot : end_frame
    % Clear figure so that it only displays figures from current datapoint
    clf;
    % Set Axes properties
    ax = axes('Units', 'normalized', 'Position', [0.01 0.16 0.98 0.76]);
    hold all;
    ax.XLim = [-500 500];
    ax.YLim = [-500 500];
    %{
    ax.TickLength = [0 0];
    ax.XTick = [0:(max_sites / 5):max_sites];
    ax.XTickLabel = {0, 0.008 * max_sites / 5, 0.008 * 2 * max_sites / 5, ...
                    0.008 * 3 * max_sites / 5, 0.008 * 4 * max_sites / 5, 0.008 * max_sites};
    ax.YTickLabel = {};
    ax.XLabel.String = 'Distance from plus-end (microns)';
    %}

    % Draw MTs
    for i_mt = 1:1:n_mts
        
        n_sites = mt_lengths(i_mt);
        %{
        mt_pos = double(filament_pos(i_mt, i_data) * site_width);
        first_pos = filament_pos(1, i_data) * site_width;
        mt_height = 8 * (i_mt - 1) * site_height;

        if (n_mts > 1)
            second_pos = filament_pos(2, i_data) * site_width + 1;
            %{
            left = first_pos;
            right = second_pos + n_sites + 1;
            center = (left + right) / 2;
            ax.XLim = [center - 65 center + 65];
            %}

            if (first_pos < second_pos)
                ax.XLim = [(first_pos - 1) (second_pos + 550 + 1)];
                %ax.XLim = [(first_pos-1) (second_pos + max_sites + 1)];
            else
                ax.XLim = [(second_pos - 1) (first_pos + 550 + 1)];
                %  ax.XLim = [(second_pos-1) (first_pos + max_sites + 1)];
            end

        else
            %ax.XLim = [first_pos-1 first_pos + 250];
            ax.XLim = [first_pos (first_pos + n_sites + 1)];
        end
        %}
       % rectangle('Position', [mt_pos mt_height (n_sites + 1) site_height], ...
        %    'FaceColor', [0.8 0.8 0.8], 'Curvature', [0 0]);
        
        plus_pos = filament_pos(:, 1, i_mt, i_data);
        minus_pos = filament_pos(:, 2, i_mt, i_data);
        
        %length = sqrt((plus_pos(1) - minus_pos(1))^2 + (plus_pos(2) - minus_pos(2))^2)
        
        line([plus_pos(1), minus_pos(1)],[plus_pos(2), minus_pos(2)], 'LineWidth', 2);
        
        for i_site = 1 : n_sites
            if(protein_ids(i_site, i_mt, i_data) ~= -1)
                disp('BOOYA')
                pos = [(i_site/n_sites - n_sites/2)*8 -25 50 50];
                 rectangle('Position',pos,'FaceColor', purple, 'Curvature', [1 1]);
            end
        end

           %{
        % Draw motors
        motor_IDs = protein_ids(:, i_mt, i_data);
        motor_head_trailing = motor_trailing(:, i_mt, i_data);

        for i_motor = 1:1:n_sites
            motor_pos = i_motor * site_width + mt_pos - 1;
            motor_height = mt_height + site_height;
            pseudo_height = motor_height + 1/2 * site_height; i
            pseudo_dx = site_width;
            mt_dx = 1;
            motor_center_y = motor_height + site_height / 2;
            motor_center_x = motor_pos + site_width / 2;
            pseudo_center_y = pseudo_height + site_height / 2;
            anchor_y = mt_height + 5 * site_height / 2;
            anchor_x = motor_pos + site_width;

            if (mod(i_mt, 2) == 0)
                motor_height = mt_height - site_height;
                pseudo_height = motor_height - 1/2 * site_height;
                mt_dx = -1;
                motor_center_y = motor_height + site_height / 2;
                pseudo_center_y = pseudo_height + site_height / 2;
                anchor_y = mt_height - 3 * site_height / 2;
            end

            if motor_IDs(i_motor) ~= -1

                if i_motor == 1

                    if motor_IDs(i_motor) == motor_IDs(i_motor + 1)
                        rectangle('Position', [motor_pos motor_height site_width site_height], ...
                            'FaceColor', blue, 'Curvature', [1 1]);
                        line([motor_center_x, anchor_x], [motor_center_y, anchor_y], ...
                            'LineWidth', 1.5, 'Color', 'black');
                    else
                        rectangle('Position', [motor_pos motor_height site_width site_height], ...
                            'FaceColor', blue, 'Curvature', [1 1]);
                        line([motor_center_x, motor_center_x], [motor_center_y, anchor_y], ...
                            'LineWidth', 1.5, 'Color', 'black');

                        if (motor_head_trailing(i_motor))
                            rectangle('Position', [motor_pos - mt_dx * pseudo_dx pseudo_height site_width site_height], ...
                                'FaceColor', blue, 'Curvature', [1 1]);
                            line([motor_center_x - mt_dx * pseudo_dx, motor_center_x], [pseudo_center_y, anchor_y], ...
                                'LineWidth', 1.5, 'Color', 'black');
                        else
                            rectangle('Position', [motor_pos + mt_dx * pseudo_dx pseudo_height site_width site_height], ...
                                'FaceColor', blue, 'Curvature', [1 1]);
                            line([motor_center_x + mt_dx * pseudo_dx, motor_center_x], [pseudo_center_y, anchor_y], ...
                                'LineWidth', 1.5, 'Color', 'black');
                        end

                    end

                elseif i_motor == max_sites

                    if motor_IDs(i_motor) == motor_IDs(i_motor - 1)
                        rectangle('Position', [motor_pos motor_height site_width site_height], ...
                            'FaceColor', blue, 'Curvature', [1 1]);
                        line([motor_center_x - site_width / 2, anchor_x - site_width / 2], ...
                            [anchor_y, motor_center_y], 'LineWidth', 1.5, 'Color', 'black');
                    else
                        rectangle('Position', [motor_pos motor_height site_width site_height], ...
                            'FaceColor', blue, 'Curvature', [1 1]);
                        line([motor_center_x, motor_center_x], [motor_center_y, anchor_y], ...
                            'LineWidth', 1.5, 'Color', 'black');

                        if (motor_head_trailing(i_motor))
                            rectangle('Position', [motor_pos - mt_dx * pseudo_dx pseudo_height site_width site_height], ...
                                'FaceColor', blue, 'Curvature', [1 1]);
                            line([motor_center_x - mt_dx * pseudo_dx, motor_center_x], [pseudo_center_y, anchor_y], ...
                                'LineWidth', 1.5, 'Color', 'black');
                        else
                            rectangle('Position', [motor_pos + mt_dx * pseudo_dx pseudo_height site_width site_height], ...
                                'FaceColor', blue, 'Curvature', [1 1]);
                            line([motor_center_x + mt_dx * pseudo_dx, motor_center_x], [pseudo_center_y, anchor_y], ...
                                'LineWidth', 1.5, 'Color', 'black');
                        end

                    end

                else

                    if motor_IDs(i_motor) == motor_IDs(i_motor + 1)
                        rectangle('Position', [motor_pos motor_height site_width site_height], ...
                            'FaceColor', blue, 'Curvature', [1 1]);
                        line([motor_center_x, anchor_x], [motor_center_y, anchor_y], ...
                            'LineWidth', 1.5, 'Color', 'black');
                    elseif motor_IDs(i_motor) == motor_IDs(i_motor - 1)
                        rectangle('Position', [motor_pos motor_height site_width site_height], ...
                            'FaceColor', blue, 'Curvature', [1 1]);
                        line([motor_center_x - site_width / 2, anchor_x - site_width / 2], ...
                            [anchor_y, motor_center_y], 'LineWidth', 1.5, 'Color', 'black');
                    else
                        rectangle('Position', [motor_pos motor_height site_width site_height], ...
                            'FaceColor', blue, 'Curvature', [1 1]);
                        line([motor_center_x, motor_center_x], [motor_center_y, anchor_y], ...
                            'LineWidth', 1.5, 'Color', 'black');

                        if (motor_head_trailing(i_motor))
                            rectangle('Position', [motor_pos - mt_dx * pseudo_dx pseudo_height site_width site_height], ...
                                'FaceColor', blue, 'Curvature', [1 1]);
                            line([motor_center_x - mt_dx * pseudo_dx, motor_center_x], [pseudo_center_y, anchor_y], ...
                                'LineWidth', 1.5, 'Color', 'black');
                        else
                            rectangle('Position', [motor_pos + mt_dx * pseudo_dx pseudo_height site_width site_height], ...
                                'FaceColor', blue, 'Curvature', [1 1]);
                            line([motor_center_x + mt_dx * pseudo_dx, motor_center_x], [pseudo_center_y, anchor_y], ...
                                'LineWidth', 1.5, 'Color', 'black');
                        end

                    end

                end

            end

        end

        % Array of xlink IDs for this MT
        xlink_IDs = xlink_data(:, i_mt, i_data);
        % Array of xlink IDs for neighbor MT
        neighb_IDs = zeros(max_sites, 1) - 1;
        neighb_mt_pos = 0;
        neighb_mt_height = 0;

        if (n_mts > 1)

            if (mod(i_mt, 2) == 0)
                neighb_IDs = xlink_data(:, i_mt - 1, i_data);
                neighb_mt_pos = filament_pos(i_mt - 1, i_data) * site_width;
                neighb_mt_height = 8 * (i_mt - 2) * site_height;
            else
                neighb_IDs = xlink_data(:, i_mt + 1, i_data);
                neighb_mt_pos = filament_pos(i_mt + 1, i_data) * site_width;
                neighb_mt_height = 8 * (i_mt) * site_height;
            end

        end

        % Scan thru MTs and plot all xlinks
        for i_xlink = 1:1:n_sites
            xlink_pos = mt_pos + i_xlink * site_width - 1;
            xlink_height = mt_height + site_height;
            xlink_center_x = xlink_pos + site_width / 2;
            xlink_center_y = xlink_height + site_height / 2;

            if (mod(i_mt, 2) == 0)
                xlink_height = mt_height - site_height;
                xlink_center_y = xlink_height - site_height / 2;
            end

            if (xlink_IDs(i_xlink) ~= -1)
                double_bound = false;

                for i_scan = -xlink_cutoff:1:xlink_cutoff
                    i_neighb = i_xlink + i_scan + mt_pos - neighb_mt_pos;

                    if i_neighb < 1
                        i_neighb = 1;
                    elseif i_neighb > n_sites
                        i_neighb = n_sites;
                    end

                    if (xlink_IDs(i_xlink) == neighb_IDs(i_neighb))
                        neighb_center_x = neighb_mt_pos + (i_neighb - 1/2) * site_width;
                        neighb_center_y = neighb_mt_height + 3 * site_height / 2;
                        neighb_height = neighb_mt_height + site_height;

                        if (mod(i_mt, 2) ~= 0)
                            neighb_center_y = neighb_mt_height + site_height / 2;
                            neighb_height = neighb_mt_height;
                            % To ensure spring is only plotted once per
                            % xlink, only plot on odd-numbered MTs
                            xa = xlink_center_x; ya = xlink_height + site_height;
                            xb = double(neighb_center_x); yb = neighb_height - site_height;
                            ne = 6; a = 6; ro = 1;
                            [xs, ys] = spring(xa, ya, xb, yb, ne, a, ro);
                            plot(xs, ys, 'LineWidth', 1, 'Color', purple);
                        end

                        double_bound = true;
                        rectangle('Position', [xlink_pos xlink_height site_width site_height], ...
                            'FaceColor', purple, 'Curvature', [0.5 0.5]);
                    end

                end

                if (double_bound == false)
                    rectangle('Position', [xlink_pos xlink_height site_width site_height], ...
                        'FaceColor', 'm', 'Curvature', [0.5 0.5]);
                    xa = xlink_center_x; ya = xlink_height + site_height;
                    xb = xlink_center_x; yb = xlink_height + 4 * site_height;

                    if (mod(i_mt, 2) == 0)
                        ya = xlink_height;
                        yb = xlink_height - 3 * site_height;
                    end

                    ne = 3; a = 3; ro = 1;
                    [xs, ys] = spring(xa, ya, xb, yb, ne, a, ro);
                    plot(xs, ys, 'LineWidth', 1, 'Color', 'm');
                end

            end

        end

        % Array of tether coords for this MT
        teth_coords = teth_data(:, i_mt, i_data);

        for i_teth = 1:1:max_sites - 1

            if (teth_coords(i_teth) ~= -1)
                end_height = (mt_height + neighb_mt_height + site_height) / 2;

                if (n_mts < 2)
                    end_height = mt_height + 5 * site_height;
                end

                start_height = mt_height + 5 * site_height / 2;

                if (mod(i_mt, 2) == 0)
                    start_height = mt_height - 3 * site_height / 2;
                end

                if (teth_coords(i_teth) ~= teth_coords(i_teth + 1))
                    start_pos = i_teth * site_width + mt_pos - 1;

                    if (i_teth > 1)

                        if (motor_IDs(i_teth) ~= motor_IDs(i_teth - 1) ...
                                && motor_IDs(i_teth) ~= motor_IDs(i_teth + 1))
                            start_pos = start_pos + site_width / 2 - 1;
                        end

                    elseif (motor_IDs(i_teth) ~= motor_IDs(i_teth + 1))
                        start_pos = start_pos + site_width / 2 - 1;
                    end

                    end_pos = teth_coords(i_teth) * site_width + (3/2) * site_width - 1;
                    xa = start_pos; ya = start_height;
                    xb = end_pos; yb = end_height;
                    ne = 10; a = 10; ro = 0.5;
                    if abs(xa - xb) <= teth_cutoff
                        [xs, ys] = spring(xa, ya, xb, yb, ne, a, ro);
                        plot(xs, ys, 'LineWidth', 1, 'Color', 'black');
                    else
                        disp(xa - xb);
                        disp(teth_coords(i_teth));
                        disp(i_teth);
                    end
                end
            end
        end
      %}
    end
    dim = [0.0105 0.62 .3 .3];
    time = (i_data - 1) * time_per_datapoint;
    %time = time - 500;
    str = sprintf('Time: %#.2f seconds', time);
    annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on');
    drawnow();
    %  frame = getframe(gcf); %(fig1); %, frame_box);
    writeVideo(v, getframe(gcf));
end

close(v);
%}