clear;
close all;
sim_name = 'run_endtag_final/endtag_final_6nM_1250_0';
%sim_name = 'endtag_shortCoop_250_0';
file_dir = '/home/shane/Projects/overlap_analysis/mgh_model';

% Open log file and parse it into param labels & their values
log_file = sprintf('%s/%s.log', file_dir, sim_name);
log = textscan(fileread(log_file),'%s %s', 'Delimiter', '=');
params = log{1,1};
values = log{1,2};
% Read in number of MTs
n_mts = str2double(values{contains(params, 'count')});
n_sites = values{contains(params, 'length')};
n_sites = sscanf(n_sites, '%i');
% Read in system params
delta_t = sscanf(values{contains(params, 'delta_t')}, '%g');
total_steps = str2double(values{contains(params, 'n_steps')});
data_threshold = sscanf(values{contains(params, 'data_threshold')}, '%g');
if any(contains(params, 'DATA_THRESHOLD') ~= 0)
   data_threshold = str2double(values{contains(params, 'DATA_THRESHOLD')});
end
n_steps = total_steps - data_threshold;
% Use max possible number of datapoints to calculate time_per_datapoint (as is done in Sim)
n_datapoints = str2double(values{contains(params, 'n_datapoints')});
time_per_datapoint = delta_t * n_steps / n_datapoints;
% Use actual recorded number of datapoints to parse thru data/etc
if any(contains(params, 'N_DATAPOINTS') ~= 0)
   n_datapoints = str2double(values{contains(params, 'N_DATAPOINTS')});
end


t_start = 0; 
dwell_time = 0.05; % in sections
dwell_points = dwell_time / time_per_datapoint; 
motor_ID = 2;
xlink_ID = 1;
tubulin_ID = 0;

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
occupany = zeros(n_sites, n_mts, n_datapoints);
occuFile = sprintf('%s/%s_occupancy.file', file_dir, sim_name);
if isfile(occuFile)
    occupancy_file = fopen(occuFile);
    occu_raw_data = fread(occupancy_file, [n_datapoints * n_mts * n_sites], '*int');
    occupancy = reshape(occu_raw_data, n_sites, n_mts, n_datapoints);
    fclose(occupancy_file);
end
% motor ID data
motor_IDs = zeros(n_sites, n_mts, n_datapoints) - 1;
motorFile = sprintf('%s/%s_motorID.file', file_dir, sim_name);
if isfile(motorFile)
    motor_data_file = fopen(motorFile);
    motor_raw_data = fread(motor_data_file, [n_mts * n_sites * n_datapoints], '*int');
    fclose(motor_data_file);
    motor_IDs = reshape(motor_raw_data, n_sites, n_mts, n_datapoints);
end
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

%% %%%%%%%%%%%25 nM motors%%%%%%%%%%%%%%%%%%
%{
% close all;
gaussAmp = 1500;
bkgLevel = 200;
noiseStd = 100;
intensityMax = gaussAmp*1.75;% + bkgLevel;

load('matrix_column_time_row_site_25_180_-4_26_1_10_14_05_00_099');

figure; imagesc(matrix); colormap gray; %axis equal

% motors in overlap
tMin = 9000;
timeUse = tMin:tMin+20;
dataMatrix = sum(matrix(timeUse,:),1);

imageMotors = imageGaussianOverlap(dataMatrix,siteLength,pixelLength,pixelPad,...
    gaussSigma,gaussAmp,bkgLevel,noiseStd,doPlot);

figure, imagesc(imageMotors);
colormap gray; axis image
set(gca,'Xtick',[]); set(gca,'Ytick',[]);
imageMotors = imageMotors/intensityMax; %convert to grayscale

%line for PRC1
lineMatrix = ones(size(dataMatrix));
gaussAmp = 2000;
imageLine = imageGaussianOverlap(lineMatrix,siteLength,pixelLength,pixelPad,...
    gaussSigma,gaussAmp,bkgLevel,noiseStd,doPlot);
imageLine = circshift(imageLine,2,1);
figure, imagesc(imageLine);
colormap gray; axis image
set(gca,'Xtick',[]); set(gca,'Ytick',[]);
imageLine = imageLine/intensityMax; %convert to grayscale

%blue channel
imageBlue = ones(size(imageLine))*bkgLevel + randn(size(imageLine))*noiseStd;
imageBlue = imageBlue/intensityMax;

%merge into RGB image
imageRGB = cat(3, imageLine, imageMotors, imageBlue);
figure; imshow(imageRGB); axis image
% set(gca,'Xtick',[]); set(gca,'Ytick',[]);
image25 = imageRGB;

% l = line([1 101],[-50 -50],'Color','k','LineWidth',4);
% set(l,'clipping','off')
% l2 = line([-10 -10],[1 201],'Color','k','LineWidth',4);
% set(l2,'clipping','off')

saveas(gca,'overlap25','png');

%% %%%%%%%%%%%100 nM motors%%%%%%%%%%%%%%%%%%
% close all;
gaussAmp = 1000;
bkgLevel = 200;
noiseStd = 100;
intensityMax = gaussAmp*2.25;% + bkgLevel;

load('matrix_column_time_row_site_100_60_-4_26_1_10_14_05_00_099');

figure; imagesc(matrix); colormap gray; %axis equal

% motors in overlap
tMin = 6000;
timeUse = tMin:tMin+20;
dataMatrix = sum(matrix(timeUse,:),1);

imageMotors = imageGaussianOverlap(dataMatrix,siteLength,pixelLength,pixelPad,...
    gaussSigma,gaussAmp,bkgLevel,noiseStd,doPlot);

figure, imagesc(imageMotors);
colormap gray; axis image
set(gca,'Xtick',[]); set(gca,'Ytick',[]);
imageMotors = imageMotors/intensityMax; %convert to grayscale

%line for PRC1
lineMatrix = ones(size(dataMatrix));
gaussAmp = 3000;
noiseStd = 100;
imageLine = imageGaussianOverlap(lineMatrix,siteLength,pixelLength,pixelPad,...
    gaussSigma,gaussAmp,bkgLevel,noiseStd,doPlot);
imageLine = circshift(imageLine,2,1);
figure, imagesc(imageLine);
colormap gray; axis image
set(gca,'Xtick',[]); set(gca,'Ytick',[]);
imageLine = imageLine/intensityMax; %convert to grayscale

%blue channel
imageBlue = ones(size(imageLine))*bkgLevel + randn(size(imageLine))*noiseStd;
imageBlue = imageBlue/intensityMax;

%merge into RGB image
imageRGB = cat(3, imageLine, imageMotors, imageBlue);
figure; imshow(imageRGB); axis image
image100 = imageRGB;

% set(gca,'Xtick',[]); set(gca,'Ytick',[]);

% l = line([1 101],[-50 -50],'Color','k','LineWidth',4);
% set(l,'clipping','off')
% l2 = line([-10 -10],[1 201],'Color','k','LineWidth',4);
% set(l2,'clipping','off')

saveas(gca,'overlap100','png');

%% put 3 images together

image25Pad = padarray(image25,[0 12/2],'replicate');
image100Pad = zeros(size(image5));
image100Pad(:,13:59,:) = image100;
image100Pad(:,1:12,:) = image100Pad(:,12:23,:);
image100Pad(:,60:66,:) = image100Pad(:,41:47,:);
image100Pad(:,66:72,:) = image100Pad(:,41:47,:);

whiteLine = repmat(ones(1,72),1,1,3);

imageMerge = cat(1,image100Pad,whiteLine,image25Pad,whiteLine,image5);
figure; imshow(imageMerge); axis image

offset = 3; 
l = line(offset+[1 34.33],offset*[1 1],'Color','w','LineWidth',4);

text(51,20,'100 nM','Color','w','FontSize',18);
text(54.5,61.5,'25 nM','Color','w','FontSize',18);
text(57.5,103,'5 nM','Color','w','FontSize',18);
saveas(gca,'overlap_final','pdf');
%}
