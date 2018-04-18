clear all

% Parameters from sim
n_datapoints = 100000;
motor_ID = 2;
mt_length = 250;
n_mts = 2;
xlink_cutoff = 7;

% File info
simName = 'presXL5';
fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
mtFileName = '%s_MTcoord.file';
motorFileName = '%s_motorID.file';
xlinkFileName = '%s_xlinkID.file';
tethFileName = '%s_tether_coord.file';
mtFile = sprintf(fileDirectory, sprintf(mtFileName, simName));
motorFile = sprintf(fileDirectory, sprintf(motorFileName, simName));
xlinkFile = sprintf(fileDirectory, sprintf(xlinkFileName, simName));
tethFile = sprintf(fileDirectory, sprintf(tethFileName, simName));

% Figure parameters (i.e., how they appear)
n_frames = 100000;
frames_per_plot = 10000;
start_frame = 001;
site_height = 1;
site_width = 1;

% Videowriter details
v = VideoWriter('newfile5.avi');
v.FrameRate = (n_frames / frames_per_plot) / 25;
open(v);
frame_box = [0 0 1545 200];

% Figure details
fig1 = figure;
set(fig1, 'Position', [0 100 1600 400])

mt_data_file = fopen(mtFile);
mt_raw_data = fread(mt_data_file, [n_mts * n_datapoints], '*double');
fclose(mt_data_file);
mt_data = reshape(mt_raw_data, n_mts, n_datapoints);

motor_data_file = fopen(motorFile);
motor_raw_data = fread(motor_data_file, [n_mts * mt_length * n_datapoints], '*int');
fclose(motor_data_file);
motor_data = reshape(motor_raw_data, mt_length, n_mts, n_datapoints);

xlink_data_file = fopen(xlinkFile);
xlink_raw_data = fread(xlink_data_file, [n_mts * mt_length * n_datapoints], '*int');
fclose(xlink_data_file);
xlink_data = reshape(xlink_raw_data, mt_length, n_mts, n_datapoints);

teth_data_file = fopen(tethFile);
teth_raw_data = fread(teth_data_file, [n_mts * mt_length * n_datapoints], '*double');
fclose(teth_data_file);
teth_data = reshape(teth_raw_data, mt_length, n_mts, n_datapoints);

% Run through all datapoints; each one is a frame in our movie
for i_data=start_frame:frames_per_plot:(start_frame + n_frames - 1)
    
    % Clear figure so that it only displays figures from current datapoint
    clf;        

    % Set Axes properties
    ax = axes('Units', 'normalized', 'Position', [0.01 0.12 0.98 0.76]);
    ax.XLim = [0 (mt_length + 1)];
    ax.YLim = [-1 11];
    ax.TickLength = [0 0];
    ax.XTick = [0:(mt_length/5):mt_length];
    ax.XTickLabel = {0, 0.008*mt_length/5, 0.008*2*mt_length/5, ... 
        0.008*3*mt_length/5, 0.008*4*mt_length/5, 0.008*mt_length};
    ax.YTickLabel = {};
    ax.XLabel.String = 'Distance across MT relative to minus-end (microns)';
    
    hold all;
    
    % Draw MTs
    for i_mt=1:1:n_mts
        
        mt_pos = mt_data(i_mt, i_data)*site_width;
        first_pos = mt_data(1, i_data)*site_width;
        second_pos = mt_data(2, i_data)*site_width;
        mt_height = 8*(i_mt - 1)*site_height;
        
        if(first_pos < second_pos)
            ax.XLim = [first_pos first_pos + (3/2)*mt_length];
        else
            ax.XLim = [second_pos second_pos + (3/2)*mt_length];
        end
        
        rectangle('Position', [mt_pos mt_height (mt_length + 1) site_height], ...
            'FaceColor', [0.8 0.8 0.8], 'Curvature', [0 0]);
        
        % Draw motors
        motor_IDs = motor_data(:, i_mt, i_data);
        for i_motor=1:1:mt_length
            motor_pos = i_motor*site_width + mt_pos;
            motor_height = mt_height + site_height;
            pseudo_height = motor_height + 2*site_height;
            motor_center_y = motor_height + site_height/2;
            motor_center_x = motor_pos + site_width/2;
            pseudo_center = mt_height + 7*site_height/2;
            anchor_y = mt_height + 5*site_height/2;
            anchor_x = motor_pos + site_width;
            if(mod(i_mt, 2) == 0)
                motor_height = mt_height - site_height;
                pseudo_height = motor_height - 2*site_height;
                motor_center_y = motor_height + site_height/2;
                pseudo_center = mt_height - 5*site_height/2;
                anchor_y = mt_height - 3*site_height/2;
            end
            if motor_IDs(i_motor) ~= -1
                if i_motor == 1
                    if motor_IDs(i_motor) == motor_IDs(i_motor + 1)
                        rectangle('Position', [motor_pos motor_height site_width site_height], ...
                            'FaceColor', 'g', 'Curvature', [1 1]);
                        line([motor_center_x, anchor_x], [motor_center_y, anchor_y],  ...
                            'LineWidth', 1.5, 'Color', 'black');
                    else
                        rectangle('Position', [motor_pos motor_height site_width site_height], ...
                            'FaceColor', 'r', 'Curvature', [1 1]);
                        rectangle('Position', [motor_pos pseudo_height site_width site_height], ...
                            'FaceColor', 'r', 'Curvature', [1 1]);
                        line([motor_center_x, motor_center_x], [motor_center_y, pseudo_center], ...
                            'LineWidth', 1.5, 'Color','black');
                    end
                elseif i_motor == mt_length
                    if motor_IDs(i_motor) == motor_IDs(i_motor - 1)
                        rectangle('Position', [motor_pos motor_height site_width site_height], ...
                            'FaceColor', 'g', 'Curvature', [1 1]);
                        line([motor_center_x - site_width/2, anchor_x - site_width/2], ...
                            [anchor_y, motor_center_y], 'LineWidth', 1.5, 'Color', 'black');
                    else
                        rectangle('Position', [motor_pos motor_height site_width site_height], ...
                            'FaceColor', 'r', 'Curvature', [1 1]);
                        rectangle('Position', [motor_pos pseudo_height site_width site_height], ...
                            'FaceColor', 'r', 'Curvature', [1 1]);
                        line([motor_center_x, motor_center_x], [motor_center_y, pseudo_center], ...
                            'LineWidth', 1.5, 'Color','black');
                    end
                else
                    if motor_IDs(i_motor) == motor_IDs(i_motor + 1)
                        rectangle('Position', [motor_pos motor_height site_width site_height], ...
                            'FaceColor', 'g', 'Curvature', [1 1]);
                        line([motor_center_x, anchor_x], [motor_center_y, anchor_y],  ...
                            'LineWidth', 1.5, 'Color', 'black');
                    elseif motor_IDs(i_motor) == motor_IDs(i_motor - 1)
                        rectangle('Position', [motor_pos motor_height site_width site_height], ...
                            'FaceColor', 'g', 'Curvature', [1 1]);
                        line([motor_center_x - site_width/2, anchor_x - site_width/2], ...
                            [anchor_y, motor_center_y], 'LineWidth', 1.5, 'Color', 'black');
                    else
                        rectangle('Position', [motor_pos motor_height site_width site_height], ...
                            'FaceColor', 'r', 'Curvature', [1 1]);
                        rectangle('Position', [motor_pos pseudo_height site_width site_height], ...
                            'FaceColor', 'r', 'Curvature', [1 1]);
                        line([motor_center_x, motor_center_x], [motor_center_y, pseudo_center], ...
                            'LineWidth', 1.5, 'Color','black');
                    end
                end
            end
        end
        
        % Array of xlink IDs for this MT
        xlink_IDs = xlink_data(:, i_mt, i_data);
        % Array of xlink IDs for neighbor MT
        neighb_IDs = zeros(mt_length);
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
        for i_xlink=1:1:mt_length
            xlink_pos = mt_pos + i_xlink*site_width;
            xlink_height = mt_height + site_height;
            xlink_center_x = xlink_pos + site_width/2;
            xlink_center_y = xlink_height + site_height/2;
            if(mod(i_mt, 2) == 0)
                xlink_height = mt_height - site_height;
                xlink_center_y = xlink_height - site_height/2;
            end
            if(xlink_IDs(i_xlink) ~= -1)
                double_bound = false;
                for i_scan = -xlink_cutoff:1:xlink_cutoff
                    i_neighb = i_xlink + i_scan + mt_pos - neighb_mt_pos;
                    if i_neighb < 1
                        i_neighb = 1;
                    elseif i_neighb > mt_length
                        i_neighb = mt_length;
                    end
                    if(xlink_IDs(i_xlink) == neighb_IDs(i_neighb))
                        neighb_center_x = neighb_mt_pos + (i_neighb + 1/2)*site_width;
                        neighb_center_y = neighb_mt_height + 3*site_height/2;
                        neighb_height = neighb_mt_height + site_height;
                        if(mod(i_mt, 2) ~= 0)
                            neighb_center_y = neighb_mt_height + site_height/2;
                            neighb_height = neighb_mt_height;
                            % To ensure spring is only plotted once per
                            % xlink, only plot on odd-numbered MTs
                            xa = xlink_center_x; ya = xlink_height + site_height;
                            xb = neighb_center_x; yb = neighb_height - site_height;
                            ne = 3; a = 8*site_height; ro = 1;
                            [xs,ys] = spring(xa,ya,xb,yb,ne,a,ro);
                            plot(xs,ys,'LineWidth',2, 'Color', 'm');
                        end
                        double_bound = true;
                        rectangle('Position', [xlink_pos xlink_height site_width site_height], ...
                            'FaceColor', 'm', 'Curvature', [0.5 0.5]);
                    end
                end
                if(double_bound == false)
                    rectangle('Position', [xlink_pos xlink_height site_width site_height], ...
                        'FaceColor', 'cyan', 'Curvature', [0.5 0.5]);
                end
            end
        end
    end
    %{
    rectangle('Position', [(i_data)*site_width site_height site_width site_height], ...
        'FaceColor', 'g', 'Curvature', [1 1]);
    rectangle('Position', [(i_data + 1)*site_width site_height site_width site_height], ...
        'FaceColor', 'g', 'Curvature', [1 1]);
    plot([i_data + site_width/2, i_data + site_width], ...
        [site_height*(3/2), site_height*(5/2)], 'LineWidth', 2);
    plot([i_data + site_width*(3/2), i_data + site_width], ...
        [site_height*(3/2), site_height*(5/2)], 'LineWidth', 2);
    xa = i_data + site_width; ya = site_height*(5/2);
    xb = i_data; yb = 5 + site_height*(5/2);
    ne = 5; a = 10; ro = 0.005;
    [xs,ys] = spring(xa,ya,xb,yb,ne,a,ro); plot(xs,ys,'LineWidth',2)
   
    %}
    
    frame = getframe(fig1); %, frame_box);
    writeVideo(v, frame);
end

close(v);