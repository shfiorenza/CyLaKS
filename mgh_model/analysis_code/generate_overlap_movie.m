clear;
close all;
simName = 'Endtag_1750';
n_sites = 1750;
n_mts = 1;
n_datapoints = 10000;
n_steps = 60000000;
delta_t = 0.000025; % seconds
movie_name = "test_mov";
movie_frames_per_plot = 50;
movie_duration = 30; % real life seconds
movie_dwell_time = 5; % in sim-seconds; time between each 'image'
movie_start = 01; % datapoint to start at
time_per_frame = n_steps * delta_t / n_datapoints;

motor_ID = 2;
xlink_ID = 1;
tubulin_ID = 0;

% parameters for making simulated image (i.e., each frame)
siteLength = 8;
pixelLength = 150;
pixelPad = 20;
gaussSigma = 1.0;
doPlot = 0;
gaussAmp = 4000;
bkgLevel = 200;
noiseStd = 100;
intensityMax = gaussAmp/2;% + bkgLevel;

fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
occu_fileStruct = '%s_occupancy.file';
mt_fileStruct = '%s_mt_coord.file';
motorID_fileStruct = '%s_motorID.file';

v = VideoWriter(movie_name);
v.FrameRate = (n_datapoints / movie_frames_per_plot) / movie_duration;
open(v);
frame_box = [0 0 1000 300];

fig1 = figure;
set(fig1, 'Position', [250 300 1000 300]);

% microtubule coordinates - specifically, left edge of MT
mt_fileName = sprintf(fileDirectory, sprintf(mt_fileStruct, simName));
mt_data_file = fopen(mt_fileName);
mt_coords = fread(mt_data_file, [n_datapoints, n_mts], '*int');
fclose(mt_data_file);
% occupancy data on each MT
occu_fileName = sprintf(fileDirectory, sprintf(occu_fileStruct, simName));
occu_data_file = fopen(occu_fileName);
occupancy_raw_data = fread(occu_data_file, (n_datapoints * n_mts * n_sites), '*int');
fclose(occu_data_file);
occupancy = reshape(occupancy_raw_data, n_sites, n_mts, n_datapoints);
% motor ID data
motor_fileName = sprintf(fileDirectory, sprintf(motorID_fileStruct, simName));
motor_data_file = fopen(motor_fileName);
motor_raw_data = fread(motor_data_file, (n_mts * n_sites * n_datapoints), '*int');
fclose(motor_data_file);
motor_IDs = reshape(motor_raw_data, n_sites, n_mts, n_datapoints);

% Convert occupancy data to binary logic array of motor occupancy
occupancy(occupancy ~= motor_ID) = 0;
occupancy(occupancy == motor_ID) = 1;
matrix = occupancy;

% Run through movie frames and create each one
for i_data=movie_start:movie_frames_per_plot:n_datapoints-movie_frames_per_plot
    % Clear figure
    clf;
    hold all
    
    % green channel - motors
    dataMatrix = sum(matrix(:, :, i_data:i_data+movie_frames_per_plot), 3)';
    imageMotors = imageGaussianOverlap(dataMatrix,siteLength,pixelLength,pixelPad,...
        gaussSigma,gaussAmp,bkgLevel,noiseStd,doPlot);  
    imageMotors = imageMotors/intensityMax; %convert to grayscale
    % red channel - microtubule
    lineMatrix = ones(size(dataMatrix));
    imageLine = imageGaussianOverlap(lineMatrix,siteLength,pixelLength,pixelPad,...
        gaussSigma,gaussAmp,bkgLevel,noiseStd,doPlot);
    imageLine = imageLine/intensityMax; %convert to grayscale
    % blue channel - noise
    imageBlue = ones(size(imageLine))*bkgLevel + randn(size(imageLine))*noiseStd;
    imageBlue = imageBlue/intensityMax;
    % merge into RGB image
    imageRGB = cat(3, imageLine, imageMotors, imageBlue);
    imagesc(imageRGB); axis image
    set(gca,'Xtick',[]); set(gca,'Ytick',[]);
    
    dim = [0.15 0.62 0.5 0.25];
    time = (i_data - 1) * time_per_frame;
    str = sprintf('Time: %#.2f seconds', time);
    annotation('textbox',dim,'String',str,'FitBoxToText','on', 'BackgroundColor', [1 1 1]);
    drawnow();
    writeVideo(v, getframe(gcf)); 
end
close(v);