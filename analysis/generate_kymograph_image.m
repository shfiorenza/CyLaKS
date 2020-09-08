clear;
close all;
simName = 'test';
n_sites = 1000;
n_datapoints = 10000;
n_mts = 2;

iptsetpref('ImshowBorder','tight');
fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
occu_fileStruct = '%s_occupancy.file';
mt_fileStruct = '%s_mt_coord.file';

% microtubule coordinates - specifically, left edge of MT
fileName = sprintf(fileDirectory, sprintf(mt_fileStruct, simName));
data_file = fopen(fileName);
mt_coords = fread(data_file, [n_datapoints, n_sites], '*int');
fclose(data_file);
% occupancy data on each MT
fileName = sprintf(fileDirectory, sprintf(occu_fileStruct, simName));
data_file = fopen(fileName);
occupancy = fread(data_file, [n_datapoints, n_sites], '*int');
fclose(data_file);

%matrix = zeros([n_datapoints n_sites]); 
%matrix = occupancy(:,1,:);
matrix = occupancy;
%for i_data=1:1:n_datapoints
%    
%end

%matrix = 

%load('kymograph_data_column_time_row_site');

figure; imagesc(matrix); colormap gray; %axis equal

% parameters for making image
siteLength = 8;
pixelLength = 150;
gaussSigma = 1.0;
gaussAmp = 1000;
bkgLevel = 50;
noiseStd = 25;
doPlot = 0;

% truncation parameters - center region
siteMin = 250;%500;%625;
siteMax = n_sites;%2000;%1875;
timeMin = 1;
timeSkip = 5;
timeMax = 10000;
dataMatrix = matrix(timeMin:timeSkip:timeMax,siteMin:siteMax);

image2D = imageGaussianKymograph(dataMatrix,siteLength,pixelLength,...
    gaussSigma,gaussAmp,bkgLevel,noiseStd,doPlot);

intThresh = 1500;
image2D(image2D>intThresh) = intThresh;

figure, imagesc(image2D);
colormap gray; 
whitebg(2,'k')
axis([-6 126 -50 2050]);
daspect([150 1000 1]);

%find pixels for 10 microns
lineLength = 10000/pixelLength;
l = line([1 1+lineLength],[-20 -20],'Color','w','LineWidth',3);
set(l,'clipping','off')

l2 = line([-3 -3],[1 201],'Color','w','LineWidth',3);
set(l2,'clipping','off')
set(gca,'Xtick',[]); set(gca,'Ytick',[]);
saveas(gca,'kymograph','png');


% % truncation parameters - boundary region
% siteMin = 1;
% siteMax = 1000;%1875;
% timeMin = 1;
% timeSkip = 20;
% timeMax = 10000;
% dataMatrix = matrix(timeMin:timeSkip:timeMax,siteMin:siteMax);
% 
% image2DEdge = imageGaussianKymograph(dataMatrix,siteLength,pixelLength,...
%     gaussSigma,gaussAmp,bkgLevel,noiseStd,doPlot);
% 
% figure, imagesc(image2DEdge);
% colormap gray;

