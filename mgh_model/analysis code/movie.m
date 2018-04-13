clear all

% Parameters from sim
n_datapoints = 100000;
motor_ID = 2;
mt_length = 250;

% File info
simName = 'sim4';
fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
mtFile = '%s_MTcoord.file';
motorFile = '%s_motorID.file';
xlinkFile = '%s_xlinkID.file';
motorFileName = sprintf(fileDirectory, sprintf(motorFile, simName));

% Figure parameters (i.e., how they appear)
n_frames = 250;
site_height = 1;
site_width = 1;

% Videowriter details
v = VideoWriter('newfile4.avi');
v.FrameRate = 3;
open(v);
frame_box = [0 0 1545 200];

% Figure details
fig1 = figure;
set(fig1, 'Position', [0 100 1600 400])



motor_data_file = fopen(motorFileName);
motor_raw_data = fread(motor_data_file, [mt_length, n_datapoints], '*int');
fclose(motor_data_file);

% Run through all datapoints; each one is a frame in our movie
for i_data=1:1:n_frames
    
    % Clear figure so that it only displays figures from current datapoint
    clf;
    % Set Axes properties
    ax = axes('Units', 'normalized', 'Position', [0.01 0.12 0.98 0.76]);
    ax.XLim = [100 (mt_length + 1)]; 
    ax.YLim = [-5 20];
    ax.TickLength = [0 0];
    ax.XTickLabel = {'0.8', '1.2', '1.6', '2.0'};
    ax.YTickLabel = {};
    ax.XLabel.String = 'Distance across MT relative to minus-end (microns)';
    
    hold all;
    
    % Draw MT
    rectangle('Position', [0 site_height/2 (mt_length + 1) site_height], ...
        'FaceColor', [0.8 0.8 0.8], 'Curvature', [1 0.2]);
    
    
    % Draw motors
    temp_mt = motor_raw_data(:, i_data);
    for i_motor=1:1:mt_length
        if temp_mt(i_motor) ~= -1
            if i_motor == 1
                if temp_mt(i_motor) == temp_mt(i_motor + 1)
                    rectangle('Position', [(i_motor*site_width) site_height site_width site_height], ...
                        'FaceColor', 'g', 'Curvature', [1 1]);
                    plot([i_motor + site_width/2, i_motor + site_width], ...
                        [site_height*(3/2), site_height*(5/2)], 'LineWidth', 2, ...
                        'Color', 'black');               
                else
                    rectangle('Position', [(i_motor*site_width) site_height ...
                        site_width site_height], 'FaceColor', 'r', 'Curvature', [1 1]);
                    rectangle('Position', [(i_motor*site_width) 3*site_height ...
                        site_width site_height], 'FaceColor', 'r', 'Curvature', [1 1]);
                    plot([i_motor + site_width/2, i_motor + site_width/2], ...
                        [(3/2)*site_height, (7/2)*site_height], 'LineWidth', 2, ...
                        'Color', 'black');
                end
            elseif i_motor == mt_length
                if temp_mt(i_motor) == temp_mt(i_motor - 1)
                    rectangle('Position', [(i_motor*site_width) site_height site_width site_height], ...
                        'FaceColor', 'g', 'Curvature', [1 1]);
                    plot([(i_motor -1) + site_width*(3/2), (i_motor - 1) + site_width], ...
                        [site_height*(3/2), site_height*(5/2)], 'LineWidth', 2, ...
                        'Color', 'black');
                else
                    rectangle('Position', [(i_motor*site_width) site_height ...
                        site_width site_height], 'FaceColor', 'r', 'Curvature', [1 1]);
                    rectangle('Position', [(i_motor*site_width) 3*site_height ...
                        site_width site_height], 'FaceColor', 'r', 'Curvature', [1 1]);
                    plot([i_motor + site_width*(3/2), i_motor + site_width/2], ...
                        [(3/2)*site_height, (7/2)*site_height], 'LineWidth', 2, ...
                        'Color', 'black');
                end
            else  
                if temp_mt(i_motor) == temp_mt(i_motor + 1)
                    rectangle('Position', [i_motor site_height site_width site_height], ...
                        'FaceColor', 'g', 'Curvature', [1 1]);
                    plot([i_motor + site_width/2, i_motor + site_width], ...
                        [site_height*(3/2), site_height*(5/2)], 'LineWidth', 2, ...
                        'Color', 'black');
                elseif temp_mt(i_motor) == temp_mt(i_motor - 1)
                    rectangle('Position', [i_motor site_height site_width site_height], ...
                        'FaceColor', 'g', 'Curvature', [1 1]);
                    plot([(i_motor -1) + site_width*(3/2), (i_motor - 1) + site_width], ...
                        [site_height*(3/2), site_height*(5/2)], 'LineWidth', 2, ...
                        'Color', 'black');
                else
                    rectangle('Position', [(i_motor*site_width) site_height ...
                        site_width site_height], 'FaceColor', 'r', 'Curvature', [1 1]);
                    rectangle('Position', [(i_motor*site_width) 3*site_height ...
                        site_width site_height], 'FaceColor', 'r', 'Curvature', [1 1]);
                    plot([i_motor + site_height/2, i_motor + site_height/2], ...
                        [(3/2)*site_height, (7/2)*site_height], 'LineWidth', 2, ...
                        'Color', 'black');
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
    
    frame = getframe(fig1, frame_box);
    writeVideo(v, frame);
end

close(v);