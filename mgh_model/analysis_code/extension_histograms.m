clear all
simName = 'test';
% Parameters from sim
n_datapoints = 10000;
motor_cutoff = 18;
motor_r_0 = 14.5;
xlink_cutoff = 5;
xlink_r_0 = 0;

n_steps = 400000000;
delta_t = 0.0000025;

% File info
fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
motorExtFileName = '%s_motor_extension.file';
xlinkExtFileName = '%s_xlink_extension.file';
motorExtFile = sprintf(fileDirectory, sprintf(motorExtFileName, simName));
xlinkExtFile = sprintf(fileDirectory, sprintf(xlinkExtFileName, simName));

motor_data_file = fopen(motorExtFile);
motor_data = fread(motor_data_file, [(2*motor_cutoff + 1), n_datapoints], 'int');
fclose(motor_data_file);
xlink_data_file = fopen(xlinkExtFile);
xlink_data = fread(xlink_data_file, [xlink_cutoff + 1, n_datapoints], 'int');
fclose(xlink_data_file);

% Figure parameters (i.e., how they appear)
n_frames = 10000;
frames_per_plot = 100;
start_frame = 0001;

% Videowriter details
v = VideoWriter('histograms.avi');
v.FrameRate = (n_frames / frames_per_plot) / 60;
open(v);
frame_box = [0 0 1545 200];

% Figure details
fig1 = figure;
set(fig1, 'Position', [0 100 1600 400])

time_per_frame = delta_t * (n_steps / n_frames);

% Run through all datapoints; each one is a frame in our movie
for i_data=start_frame:frames_per_plot:n_frames
    
    % Clear figure so that it only displays figures from current datapoint
    clf;      
    
    subplot(2, 1, 1)
    bar(motor_data(:, i_data), 'XData', 0:1:2*motor_cutoff);
    line([2*motor_r_0 2*motor_r_0], [0 50], 'LineStyle', '--', 'Color', 'red');
    title('Extension histogram for motors');
   % ylim([0 20]);
    ylabel('Occupancy');
    xlabel('Horziontal distance doubled (sites)');
   
    subplot(2, 1, 2)
    bar(xlink_data(:, i_data), 'XData', 0:1:xlink_cutoff)
    line([xlink_r_0 xlink_r_0], [0 50], 'LineStyle', '--', 'Color', 'red');
    title('Extension histogram for crosslinks');
   % ylim([0 50]);
    ylabel('Occupancy');
    xlabel('Horizontal distance (sites)');
        
    dim = [0.8 0.23 .3 .3];
    time = (i_data - 1) * time_per_frame;
    str = sprintf('Time: %#.3g seconds', time);
    annotation('textbox',dim,'String',str,'FitBoxToText','on');
    
    frame = getframe(fig1);
    writeVideo(v, frame);
end

close(v);