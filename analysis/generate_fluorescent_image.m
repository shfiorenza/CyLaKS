clear variables;
close all;

sim_name = 'test';
%sim_name = 'run_endtag_final/endtag_final_6nM_1250_0';

t_start = 0; 
dwell_time = 0.05; % in sections

% parameters for making image
siteLength = 8;
pixelLength = 150;
pixelPad = 10;
gaussSigma = 2.0;
doPlot = 0;
gaussAmp = 1000;
bkgLevel = 20;
noiseStd = 50;
intensityMax = gaussAmp/2  + bkgLevel;

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

dwell_points = dwell_time / time_per_datapoint; 

tubulin_ID = 0;
xlink_ID = 1;
motor_ID = 2;

% microtubule coordinates - specifically, left edge of MT
mt_data = zeros(n_mts, n_datapoints);
mtFile = sprintf('%s/%s_mt_coord.file', fileDirectory, sim_name);
if isfile(mtFile)
    mt_data_file = fopen(mtFile);
    mt_raw_data = fread(mt_data_file, [n_mts * n_datapoints], '*int');
    fclose(mt_data_file);
    mt_data = reshape(mt_raw_data, n_mts, n_datapoints);
end
% occupancy data on each MT
occupancy = zeros(max_sites, n_mts, n_datapoints);
occuFile = sprintf('%s/%s_occupancy.file', fileDirectory, sim_name);
if isfile(occuFile)
    occupancy_file = fopen(occuFile);
    occu_raw_data = fread(occupancy_file, [n_datapoints * n_mts * max_sites], '*int');
    occupancy = reshape(occu_raw_data, max_sites, n_mts, n_datapoints);
    fclose(occupancy_file);
end
% motor ID data
motor_IDs = zeros(max_sites, n_mts, n_datapoints) - 1;
motorFile = sprintf('%s/%s_motorID.file', fileDirectory, sim_name);
if isfile(motorFile)
    motor_data_file = fopen(motorFile);
    motor_raw_data = fread(motor_data_file, [n_mts * max_sites * n_datapoints], '*int');
    fclose(motor_data_file);
    motor_IDs = reshape(motor_raw_data, max_sites, n_mts, n_datapoints);
end

% Plot motors
occupancy(occupancy ~= motor_ID) = 0;
occupancy(occupancy == motor_ID) = 1;
matrix = double(occupancy);
dataMatrix = sum(matrix(:, :, (n_datapoints - dwell_points - t_start/time_per_datapoint): n_datapoints - t_start/time_per_datapoint), 3)';
imageMotors = imageGaussianOverlap(dataMatrix,siteLength,pixelLength,pixelPad,...
    gaussSigma,gaussAmp,bkgLevel,noiseStd,doPlot);
imageMotors = imageMotors/intensityMax; %convert to grayscale

% line for microtubule ??
lineMatrix = ones(size(dataMatrix));
gaussAmp = 4000;
noiseStd = 125;
intensityMax = gaussAmp/2  + bkgLevel;
imageLine = imageGaussianOverlap(lineMatrix,siteLength,pixelLength,pixelPad,...
    gaussSigma,gaussAmp,bkgLevel,noiseStd,doPlot);
imageLine = imageLine/intensityMax; %convert to grayscale
% blue channel
imageBlue = ones(size(imageLine))*bkgLevel + randn(size(imageLine))*noiseStd;
imageBlue = imageBlue/intensityMax;

% merge into RGB image
imageRGB = cat(3, imageLine, imageMotors, imageBlue);
figure; imagesc(imageRGB); axis image
set(gca,'Xtick',[]); set(gca,'Ytick',[]);
image5 = imageRGB;

%saveas(gca,'overlap5','png');
