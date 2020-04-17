clear all;
% Parameters from sim
mt_lengths = [11, 11]; %, 500];
max_sites = max(mt_lengths);
n_mts = length(mt_lengths);
simName = "xlink_diffusion_doubleC_0_0";
simName = "test_1_0.75";
dur_sec = 60;
% Pseudo-constant variables
n_steps = 4000000;
n_datapoints = 10000;
start_frame = 1;
end_frame = n_datapoints;
frames_per_plot = 100;
delta_t = 0.000025; 
xlink_cutoff = 5;
cutoff = 19;
% Colors
blue = [30 144 255] / 255;
purple = [128 0 128] / 255;

% File info
movie_name = 'test';
%fileDirectory = '/home/shane/Desktop/slide_scan/%s';
fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
mtFileName = '%s_mt_coord.file';
motorFileName = '%s_motorID.file';
motorHeadFileName = '%s_motor_head_status.file';
xlinkFileName = '%s_xlinkID.file';
tethFileName = '%s_tether_coord.file';
mtFile = sprintf(fileDirectory, sprintf(mtFileName, simName));
motorFile = sprintf(fileDirectory, sprintf(motorFileName, simName));
motorHeadFile = sprintf(fileDirectory, sprintf(motorHeadFileName, simName));
xlinkFile = sprintf(fileDirectory, sprintf(xlinkFileName, simName));
tethFile = sprintf(fileDirectory, sprintf(tethFileName, simName));

% Figure parameters (i.e., how they appear)
n_frames = n_datapoints;
site_height = 1;
site_width = 1;
%{
end_frame = start_frame + n_frames - 1;
if(end_frame > n_datapoints)
    end_frame = n_datapoints;
end
%}
active_frames = end_frame - start_frame;
time_per_frame = delta_t * (n_steps / n_frames);

% Videowriter details
v = VideoWriter(movie_name);
v.FrameRate = (active_frames / frames_per_plot) / dur_sec;
open(v);
frame_box = [0 0 1545 200];

% Figure details
fig1 = figure;
%set(fig1, 'Position', [0 100 1000 250])
set(fig1, 'Position', [0 100 1600 400]);

mt_data_file = fopen(mtFile);
mt_raw_data = fread(mt_data_file, [n_mts * n_datapoints], '*double');
fclose(mt_data_file);
mt_data = reshape(mt_raw_data, n_mts, n_datapoints);

motor_data_file = fopen(motorFile);
motor_raw_data = fread(motor_data_file, [n_mts * max_sites * n_datapoints], '*int');
fclose(motor_data_file);
motor_data = reshape(motor_raw_data, max_sites, n_mts, n_datapoints);

motor_head_status_file = fopen(motorHeadFile);
motor_head_raw_data = fread(motor_head_status_file, [n_mts * max_sites * n_datapoints], '*bool');
fclose(motor_head_status_file);
motor_head_data = reshape(motor_head_raw_data, max_sites, n_mts, n_datapoints);

xlink_data_file = fopen(xlinkFile);
xlink_raw_data = fread(xlink_data_file, [n_mts * max_sites * n_datapoints], '*int');
fclose(xlink_data_file);
xlink_data = reshape(xlink_raw_data, max_sites, n_mts, n_datapoints);

teth_data_file = fopen(tethFile);
teth_raw_data = fread(teth_data_file, [n_mts * max_sites * n_datapoints], '*double');
fclose(teth_data_file);
teth_data = reshape(teth_raw_data, max_sites, n_mts, n_datapoints);

% Run through all datapoints; each one is a frame in our movie
for i_data=start_frame:frames_per_plot:end_frame
    
    % Clear figure so that it only displays figures from current datapoint
    clf;        
    % Set Axes properties
    ax = axes('Units', 'normalized', 'Position', [0.01 0.16 0.98 0.76]);
    hold all;
    ax.XLim = [0 (max_sites + 1)];
    ax.YLim = [-1 11];
    ax.TickLength = [0 0];
    ax.XTick = [0:(max_sites/5):max_sites];
    ax.XTickLabel = {0, 0.008*max_sites/5, 0.008*2*max_sites/5, ... 
        0.008*3*max_sites/5, 0.008*4*max_sites/5, 0.008*max_sites};
    ax.YTickLabel = {};
    ax.XLabel.String = 'Distance from plus-end (microns)';
    
    % Draw MTs
    for i_mt=1:1:n_mts
        n_sites = mt_lengths(i_mt);
        
        mt_pos = double(mt_data(i_mt, i_data)*site_width);
        first_pos = mt_data(1, i_data)*site_width;
        mt_height = 8*(i_mt - 1)*site_height;
        if(n_mts > 1)
            second_pos = mt_data(2, i_data)*site_width + 1;
            %{
            left = first_pos;
            right = second_pos + n_sites + 1;
            center = (left + right) / 2;
            ax.XLim = [center-65 center+65];
            %}
            
            if(first_pos < second_pos)
                ax.XLim = [(first_pos-1) (second_pos + 550 + 1)];
                %ax.XLim = [(first_pos-1) (second_pos + max_sites + 1)];
            else
                ax.XLim = [(second_pos-1) (first_pos + 550 + 1)];
              %  ax.XLim = [(second_pos-1) (first_pos + max_sites + 1)];
            end
              
        else
            %ax.XLim = [first_pos-1 first_pos + 250];
            ax.XLim = [first_pos (first_pos + n_sites + 1)];
        end
        
        rectangle('Position', [mt_pos mt_height (n_sites + 1) site_height], ...
            'FaceColor', [0.8 0.8 0.8], 'Curvature', [0 0]);
        
        % Draw motors
        motor_IDs = motor_data(:, i_mt, i_data);
        motor_head_trailing = motor_head_data(:, i_mt, i_data);
        for i_motor=1:1:n_sites
            motor_pos = i_motor*site_width + mt_pos - 1;
            motor_height = mt_height + site_height;
            pseudo_height = motor_height + 1/2*site_height;
            pseudo_dx = site_width;
            mt_dx = 1;
            motor_center_y = motor_height + site_height/2;
            motor_center_x = motor_pos + site_width/2;
            pseudo_center_y = pseudo_height + site_height/2;
            anchor_y = mt_height + 5*site_height/2;
            anchor_x = motor_pos + site_width;
            if(mod(i_mt, 2) == 0)
                motor_height = mt_height - site_height;
                pseudo_height = motor_height - 1/2*site_height;
                mt_dx = -1;
                motor_center_y = motor_height + site_height/2;
                pseudo_center_y = pseudo_height + site_height/2;
                anchor_y = mt_height - 3*site_height/2;
            end
            if motor_IDs(i_motor) ~= -1
                if i_motor == 1
                    if motor_IDs(i_motor) == motor_IDs(i_motor + 1)
                        rectangle('Position', [motor_pos motor_height site_width site_height], ...
                            'FaceColor', blue, 'Curvature', [1 1]);
                        line([motor_center_x, anchor_x], [motor_center_y, anchor_y],  ...
                            'LineWidth', 1.5, 'Color', 'black');
                    else
                        rectangle('Position', [motor_pos motor_height site_width site_height], ...
                            'FaceColor', blue, 'Curvature', [1 1]);
                        line([motor_center_x, motor_center_x], [motor_center_y, anchor_y], ...
                            'LineWidth', 1.5, 'Color','black');
                        if(motor_head_trailing(i_motor))
                            rectangle('Position', [motor_pos-mt_dx*pseudo_dx pseudo_height site_width site_height], ...
                                'FaceColor', blue, 'Curvature', [1 1]);
                             line([motor_center_x-mt_dx*pseudo_dx, motor_center_x], [pseudo_center_y, anchor_y], ...
                                'LineWidth', 1.5, 'Color', 'black');
                        else
                             rectangle('Position', [motor_pos+mt_dx*pseudo_dx pseudo_height site_width site_height], ...
                                'FaceColor', blue, 'Curvature', [1 1]);
                             line([motor_center_x+mt_dx*pseudo_dx, motor_center_x], [pseudo_center_y, anchor_y], ...
                                'LineWidth', 1.5, 'Color', 'black');
                        end
                    end
                elseif i_motor == max_sites
                    if motor_IDs(i_motor) == motor_IDs(i_motor - 1)
                        rectangle('Position', [motor_pos motor_height site_width site_height], ...
                            'FaceColor', blue, 'Curvature', [1 1]);
                        line([motor_center_x - site_width/2, anchor_x - site_width/2], ...
                            [anchor_y, motor_center_y], 'LineWidth', 1.5, 'Color', 'black');
                    else
                        rectangle('Position', [motor_pos motor_height site_width site_height], ...
                            'FaceColor', blue, 'Curvature', [1 1]);
                        line([motor_center_x, motor_center_x], [motor_center_y, anchor_y], ...
                            'LineWidth', 1.5, 'Color','black');
                        if(motor_head_trailing(i_motor))
                            rectangle('Position', [motor_pos-mt_dx*pseudo_dx pseudo_height site_width site_height], ...
                                'FaceColor', blue, 'Curvature', [1 1]);
                             line([motor_center_x-mt_dx*pseudo_dx, motor_center_x], [pseudo_center_y, anchor_y], ...
                                'LineWidth', 1.5, 'Color', 'black');
                        else
                             rectangle('Position', [motor_pos+mt_dx*pseudo_dx pseudo_height site_width site_height], ...
                                'FaceColor', blue, 'Curvature', [1 1]);
                             line([motor_center_x+mt_dx*pseudo_dx, motor_center_x], [pseudo_center_y, anchor_y], ...
                                'LineWidth', 1.5, 'Color', 'black');
                        end
                    end
                else
                    if motor_IDs(i_motor) == motor_IDs(i_motor + 1)
                        rectangle('Position', [motor_pos motor_height site_width site_height], ...
                            'FaceColor', blue, 'Curvature', [1 1]);
                        line([motor_center_x, anchor_x], [motor_center_y, anchor_y],  ...
                            'LineWidth', 1.5, 'Color', 'black');
                    elseif motor_IDs(i_motor) == motor_IDs(i_motor - 1)
                        rectangle('Position', [motor_pos motor_height site_width site_height], ...
                            'FaceColor', blue, 'Curvature', [1 1]);
                        line([motor_center_x - site_width/2, anchor_x - site_width/2], ...
                            [anchor_y, motor_center_y], 'LineWidth', 1.5, 'Color', 'black');
                    else
                        rectangle('Position', [motor_pos motor_height site_width site_height], ...
                            'FaceColor', blue, 'Curvature', [1 1]);
                        line([motor_center_x, motor_center_x], [motor_center_y, anchor_y], ...
                            'LineWidth', 1.5, 'Color','black');
                        if(motor_head_trailing(i_motor))
                            rectangle('Position', [motor_pos-mt_dx*pseudo_dx pseudo_height site_width site_height], ...
                                'FaceColor', blue, 'Curvature', [1 1]);
                             line([motor_center_x-mt_dx*pseudo_dx, motor_center_x], [pseudo_center_y, anchor_y], ...
                                'LineWidth', 1.5, 'Color', 'black');
                        else
                             rectangle('Position', [motor_pos+mt_dx*pseudo_dx pseudo_height site_width site_height], ...
                                'FaceColor', blue, 'Curvature', [1 1]);
                             line([motor_center_x+mt_dx*pseudo_dx, motor_center_x], [pseudo_center_y, anchor_y], ...
                                'LineWidth', 1.5, 'Color', 'black');
                        end
                    end
                end
            end
        end
        % Array of xlink IDs for this MT
        xlink_IDs = xlink_data(:, i_mt, i_data);
        % Array of xlink IDs for neighbor MT
        neighb_IDs = zeros(max_sites);
        neighb_IDs(:) = -1;
        neighb_mt_pos = 0;
        neighb_mt_height = 0;
        if(n_mts > 1)
            if(mod(i_mt, 2) == 0)
                neighb_IDs = xlink_data(:, i_mt - 1, i_data);
                neighb_mt_pos = mt_data(i_mt - 1, i_data)*site_width;
                neighb_mt_height = 8*(i_mt - 2)*site_height;
            else
                neighb_IDs = xlink_data(:, i_mt + 1, i_data);
                neighb_mt_pos = mt_data(i_mt + 1, i_data)*site_width;
                neighb_mt_height = 8*(i_mt)*site_height;
            end
        end
        % Scan thru MTs and plot all xlinks
        for i_xlink=1:1:n_sites
            xlink_pos = mt_pos + i_xlink*site_width - 1;
            xlink_height = mt_height + site_height;
            xlink_center_x = xlink_pos + site_width/2;
            xlink_center_y = xlink_height + site_height/2;
            if(mod(i_mt, 2) == 0)
                xlink_height = mt_height - site_height;
                xlink_center_y = xlink_height - site_height/2;
            end
            if(xlink_IDs(i_xlink) ~= -1)
                %disp(i_xlink)
                double_bound = false;
                for i_scan = -(xlink_cutoff + 1) : xlink_cutoff + 1
                    i_neighb = round(i_xlink + i_scan + mt_pos - neighb_mt_pos);
                    if i_neighb < 1
                        i_neighb = 1;
                    elseif i_neighb > n_sites
                        i_neighb = n_sites;
                    end
                    if(xlink_IDs(i_xlink) == neighb_IDs(i_neighb))
                        neighb_center_x = neighb_mt_pos + (i_neighb - 1/2)*site_width;
                        neighb_center_y = neighb_mt_height + 3*site_height/2;
                        neighb_height = neighb_mt_height + site_height;
                        if(mod(i_mt, 2) ~= 0)
                            neighb_center_y = neighb_mt_height + site_height/2;
                            neighb_height = neighb_mt_height;
                            % To ensure spring is only plotted once per
                            % xlink, only plot on odd-numbered MTs
                            xa = xlink_center_x; ya = xlink_height + site_height;
                            xb = double(neighb_center_x); yb = neighb_height - site_height;
                            ne = 6; a = 6; ro = 1;
                            [xs,ys] = spring(xa,ya,xb,yb,ne,a,ro);
                            plot(xs,ys,'LineWidth', 1, 'Color', purple);
                        end
                        double_bound = true;
                        rectangle('Position', [xlink_pos xlink_height site_width site_height], ...
                            'FaceColor', purple, 'Curvature', [0.5 0.5]);
                    end
                end
                if(double_bound == false)
                    rectangle('Position', [xlink_pos xlink_height site_width site_height], ...
                        'FaceColor', 'm', 'Curvature', [0.5 0.5]);
                    xa = xlink_center_x; ya = xlink_height + site_height;
                    xb = xlink_center_x; yb = xlink_height + 4*site_height;
                    if(mod(i_mt, 2) == 0)
                       ya = xlink_height;
                       yb = xlink_height - 3*site_height; 
                    end
                    ne = 3; a = 3; ro = 1;
                    [xs,ys] = spring(xa,ya,xb,yb,ne,a,ro);
                    plot(xs,ys,'LineWidth', 1, 'Color', 'm');
                end
            end
        end    
        % Array of tether coords for this MT
        teth_coords = teth_data(:, i_mt, i_data);
        for i_teth=1:1:max_sites - 1
            if(teth_coords(i_teth) ~= -1)
                end_height = (mt_height + neighb_mt_height + site_height) / 2;
                if(n_mts < 2)
                    end_height = mt_height + 5*site_height; 
                end
                start_height = mt_height + 5*site_height / 2;
                if(mod(i_mt, 2) == 0)
                    start_height = mt_height - 3*site_height / 2;
                end
                if(teth_coords(i_teth) ~= teth_coords(i_teth + 1))
                    start_pos = i_teth*site_width + mt_pos - 1;
                    if(i_teth > 1)
                        if(motor_IDs(i_teth) ~= motor_IDs(i_teth - 1) ...
                        && motor_IDs(i_teth) ~= motor_IDs(i_teth + 1))
                            start_pos = start_pos + site_width/2 - 1;
                        end
                    elseif(motor_IDs(i_teth) ~= motor_IDs(i_teth + 1))
                        start_pos = start_pos + site_width/2 - 1;
                    end                    
                    end_pos = teth_coords(i_teth)*site_width + (3/2)*site_width - 1;
                    xa = start_pos; ya = start_height;
                    xb = end_pos; yb = end_height;
                    ne = 10; a = 10; ro = 0.5;
                    if abs(xa - xb) <= cutoff
                        [xs,ys] = spring(xa,ya,xb,yb,ne,a,ro);
                        plot(xs,ys,'LineWidth', 1, 'Color', 'black');
                    else
                        disp(xa - xb);
                        disp(teth_coords(i_teth));
                        disp(i_teth);
                    end
               end
            end 
        end
    end
    dim = [0.0105 0.62 .3 .3];
    time = (i_data - 1) * time_per_frame;
    %time = time - 500;
    str = sprintf('Time: %#.2f seconds', time);
    annotation('textbox',dim,'String',str,'FitBoxToText','on');
    drawnow();
  %  frame = getframe(gcf); %(fig1); %, frame_box);
    writeVideo(v, getframe(gcf));
end

close(v);