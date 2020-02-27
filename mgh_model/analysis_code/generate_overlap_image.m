clear;
close all;
simName = 'EndtagC_500';
n_mts = 1;
n_sites = 500;
delta_t = 0.000025;
n_steps = 60000000;
n_datapoints = 10000;
time_per_datapoint = delta_t * n_steps / n_datapoints;
dwell_time = 6; % in sections
dwell_points = dwell_time / time_per_datapoint; 

motor_ID = 2;
xlink_ID = 1;
tubulin_ID = 0;

fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
occu_fileStruct = '%s_occupancy.file';
mt_fileStruct = '%s_mt_coord.file';
motorID_fileStruct = '%s_motorID.file';

% microtubule coordinates - specifically, left edge of MT
mt_fileName = sprintf(fileDirectory, sprintf(mt_fileStruct, simName));
mt_data_file = fopen(mt_fileName);
mt_coords = fread(mt_data_file, [n_datapoints, n_mts], '*int');
fclose(mt_data_file);
% occupancy data on each MT
occu_fileName = sprintf(fileDirectory, sprintf(occu_fileStruct, simName));
occu_data_file = fopen(occu_fileName);
occupancy_raw_data = fread(occu_data_file, [n_datapoints * n_mts * n_sites], '*int');
fclose(occu_data_file);
occupancy = reshape(occupancy_raw_data, n_sites, n_mts, n_datapoints);
% motor ID data
motor_fileName = sprintf(fileDirectory, sprintf(motorID_fileStruct, simName));
motor_data_file = fopen(motor_fileName);
motor_raw_data = fread(motor_data_file, [n_mts * n_sites * n_datapoints], '*int');
fclose(motor_data_file);
motor_IDs = reshape(motor_raw_data, n_sites, n_mts, n_datapoints);


% parameters for making image
siteLength = 8;
pixelLength = 150;
pixelPad = 20;
gaussSigma = 1.0;
doPlot = 0;

gaussAmp = 4000;
bkgLevel = 200;
noiseStd = 100;
intensityMax = gaussAmp/2;% + bkgLevel;

% Plot motors

occupancy(occupancy ~= motor_ID) = 0;
occupancy(occupancy == motor_ID) = 1;
matrix = occupancy;

% motors in overlap
dataMatrix = sum(matrix(:, :, (n_datapoints - dwell_points):n_datapoints), 3)';

imageMotors = imageGaussianOverlap(dataMatrix,siteLength,pixelLength,pixelPad,...
    gaussSigma,gaussAmp,bkgLevel,noiseStd,doPlot);

%{
figure, imagesc(imageMotors);
colormap gray; axis image
set(gca,'Xtick',[]); set(gca,'Ytick',[]);
%}
imageMotors = imageMotors/intensityMax; %convert to grayscale

% line for microtubule ??
lineMatrix = ones(size(dataMatrix));
%gaussAmp = 5000;
%noiseStd = 125;
imageLine = imageGaussianOverlap(lineMatrix,siteLength,pixelLength,pixelPad,...
    gaussSigma,gaussAmp,bkgLevel,noiseStd,doPlot);
%imageLine = circshift(imageLine,2,1);
%{
figure, hold all
imagesc(imageLine);
colormap gray; axis image
set(gca,'Xtick',[]); set(gca,'Ytick',[]);
%}
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
