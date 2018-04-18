clear all

% Parameters from sim
n_datapoints = 100000;
motor_ID = 2;
mt_length = 1250;
n_mts = 2;

% File info
simName = 'test';
fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
mtFile = '%s_MTcoord.file';
motorFile = '%s_motorID.file';
xlinkFile = '%s_xlinkID.file';
mtFileName = sprintf(fileDirectory, sprintf(mtFile, simName));
motorFileName = sprintf(fileDirectory, sprintf(motorFile, simName));

% Figure parameters (i.e., how they appear)
n_frames = 1000;
site_height = 1;
site_width = 1;

% Videowriter details
v = VideoWriter('newfile5.avi');
v.FrameRate = (n_frames / 25);
open(v);
frame_box = [0 0 1545 200];

% Figure details
fig1 = figure;
set(fig1, 'Position', [0 100 1600 500])

mt_data_file = fopen(mtFileName);
mt_raw_data = fread(mt_data_file, [n_mts, n_datapoints], '*int');
fclose(mt_data_file);

for i_data = 1:1:n_frames
   
    % Clear figure so that it only displays figures from current datapoint
    clf;
    % Set Axes properties
    ax = axes('Units', 'normalized', 'Position', [0.01 0.12 0.98 0.76]);
    ax.YLim = [0 4];
    ax.XLim = [(-mt_length/5) (2*mt_length)];
    
    temp_coords = mt_raw_data(:, i_data);
    
    for i_mt = 1:1:n_mts
        
        mt_pos = temp_coords(i_mt);
        
        rectangle('Position', [mt_pos (site_height/2 + i_mt) (mt_length + 1) site_height/2], ...
            'FaceColor', [0.8 0.8 0.8], 'Curvature', [1 0.2]);
    
    end
    
    frame = getframe(fig1);
    writeVideo(v, frame);
end

close(v);



