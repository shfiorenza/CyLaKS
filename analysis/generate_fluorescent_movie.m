clear;
close all;

% FIXME: this script needs MT positions to be updated
% Old: coordinates are scalar; x-pos of left-most edge of MTs
% New: coordinates are vectors; x,y-pos of MT center of mass 

sim_name = 'run_mobility_both/mobility_both_80_0';
movie_name = "movie_final_80pM_end";

i_start = 45000; % datapoint to start at
i_end = 50000;
datapoints_per_plot = 10; % How often to plot the fluorscent image
datapoints_dwell_time = 1; % Dwell time of virtual 'camera' creating image
movie_duration = 60; % real life seconds
scale_factor = 10;  % Controls how bright a single motor is (1 typically)

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

% Open log file and parse it into param labels & their values
fileDirectory = '../%s';
log_file = sprintf(fileDirectory, sprintf('%s.log', sim_name));
log = textscan(fileread(log_file), '%s %s', 'Delimiter', '=');
params = log{1, 1};
values = log{1, 2};
% Read in system params
dt = sscanf(values{contains(params, "dt ")}, '%g');
steps_per_datapoint = str2double(values{contains(params, "n_steps_per_snapshot ")});
time_per_datapoint = dt * steps_per_datapoint;
n_datapoints = str2double(values{contains(params, "n_datapoints ")});
% Use actual recorded number of datapoints to parse thru data/etc
if any(contains(params, "N_DATAPOINTS ") ~= 0)
    n_datapoints = str2double(values{contains(params, "N_DATAPOINTS ")});
end
site_size =  sscanf(values{contains(params, "site_size ")}, '%g') / 1000; % in um
% Read in number of MTs
n_mts = sscanf(values{contains(params, "count ")}, '%g');
if any(contains(params, "COUNT ") ~= 0)
    n_mts = sscanf(values{contains(params, "COUNT ")}, '%g');
end
% Read in MT lengths (in n_sites)
mt_lengths = zeros(1, n_mts);
for i_mt = 1 : n_mts
    string = sprintf("n_sites[%i] ", i_mt - 1);
    mt_lengths(i_mt) = sscanf(values{contains(params, string)}, '%i');
    if any(contains(params, sprintf("N_SITES[%i] ", i_mt - 1)) ~= 0)
        string = sprintf("N_SITES[%i] ", i_mt - 1);
        mt_lengths(i_mt) = sscanf(values{contains(params, string)}, '%i');
    end
end
max_sites = max(mt_lengths);

v = VideoWriter(movie_name);
v.FrameRate = ((i_end - i_start) / datapoints_per_plot) / movie_duration;
open(v);
frame_box = [0 0 1000 300];

fig1 = figure;
set(fig1, 'Position', [250 300 1000 300]);

% microtubule coordinates - specifically, left edge of MT
mt_data = zeros(n_mts, n_datapoints);
mtFile = sprintf('%s/%s_mt_coord.file', file_dir, sim_name);
if isfile(mtFile)
    mt_data_file = fopen(mtFile);
    mt_raw_data = fread(mt_data_file, [n_mts * n_datapoints], '*int');
    fclose(mt_data_file);
    mt_data = reshape(mt_raw_data, n_mts, n_datapoints);
end
% occupancy data on each MT
occupany = zeros(max_sites, n_mts, n_datapoints);
occuFile = sprintf('%s/%s_occupancy.file', file_dir, sim_name);
if isfile(occuFile)
    occupancy_file = fopen(occuFile);
    occu_raw_data = fread(occupancy_file, [n_datapoints * n_mts * max_sites], '*int');
    occupancy = reshape(occu_raw_data, max_sites, n_mts, n_datapoints);
    fclose(occupancy_file);
end
% motor ID data
motor_IDs = zeros(max_sites, n_mts, n_datapoints) - 1;
motorFile = sprintf('%s/%s_motorID.file', file_dir, sim_name);
if isfile(motorFile)
    motor_data_file = fopen(motorFile);
    motor_raw_data = fread(motor_data_file, [n_mts * max_sites * n_datapoints], '*int');
    fclose(motor_data_file);
    motor_IDs = reshape(motor_raw_data, max_sites, n_mts, n_datapoints);
end

tubulin_ID = 0;
xlink_ID = 1;
motor_ID = 2;
% Convert occupancy data to binary logic array of motor occupancy
occupancy(occupancy ~= motor_ID) = 0;
occupancy(occupancy == motor_ID) = 1;
matrix = scale_factor * occupancy;

% Run through movie frames and create each one
for i_data = i_start : datapoints_per_plot: (i_end - datapoints_dwell_time)
    % Clear figure
    clf;
    hold all
    
    % green channel - motors
    dataMatrix = sum(matrix(:, :, i_data:i_data + datapoints_dwell_time), 3)';
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
    
    dim = [0.21 0.62 0.4 0.25];
    time = (i_data - 1) * time_per_datapoint;
    str = sprintf('Time: %i seconds', int32(time));
    annotation('textbox',dim,'String',str,'FitBoxToText','on', 'BackgroundColor', [1 1 1]);
    drawnow();
    writeVideo(v, getframe(gcf)); 
end
close(v);