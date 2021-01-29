
clear variables; 

file_dir = '/home/shane/projects/CyLaKS';
sim_name = 'run_hybrid_motor/hybrid_motor_0.05_0';
%sim_name = 'ablation';
%sim_name = 'run_hetero_tubulin/hetero_tubulin_1_0';
%sim_name = 'run_endtag_vs_coop/endtag_1750_1000_0';
dwell_time = 0.1;  % dwell time of theoretical camera

i_start = 4150; %1150; %4150;
i_end = 4650; %1650; %4650;

scale_x = 10; % microns
scale_t = 5; % seconds

% parameters for making simulated image (i.e., each frame)
scale_factor = 50; % 15 % Controls how bright a single motor is (1 typically)
siteLength = 8;
pixelLength = 150; % 60; %150; % 15;
pixelPad = 5;
gaussSigma = 1.0;
doPlot = 0;
gaussAmp = 4000;
bkgLevel = 200;
noiseStd = 100;
intensityMax = gaussAmp/2 + bkgLevel;

% Open log file and parse it into param labels & their values
log_file = sprintf('%s/%s.log', file_dir, sim_name);
log = textscan(fileread(log_file), '%s %s', 'Delimiter', '=');
params = log{1, 1};
values = log{1, 2};
% Read in number of MTs
n_mts = sscanf(values{contains(params, "count ")}, '%g');
if any(contains(params, "COUNT ") ~= 0)
    n_mts = sscanf(values{contains(params, "COUNT ")}, '%g');
end
mt_lengths = zeros(1, n_mts);
for i_mt = 1 : n_mts
    string = sprintf("n_sites[%i] ", i_mt - 1);
    mt_lengths(i_mt) = sscanf(values{contains(params, string)}, '%i');
    if any(contains(params, sprintf("N_SITES[%i] ", i_mt - 1)) ~= 0)
        string = sprintf("N_SITES[%i] ", i_mt - 1);
        mt_lengths(i_mt) = sscanf(values{contains(params, string)}, '%i');
    end
end
% Read in system params
dt = sscanf(values{contains(params, "dt ")}, '%g');
steps_per_datapoint = str2double(values{contains(params, "n_steps_per_snapshot ")});
time_per_datapoint = dt * steps_per_datapoint;
n_datapoints = str2double(values{contains(params, "n_datapoints ")});
% Use actual recorded number of datapoints to parse thru data/etc
if any(contains(params, "N_DATAPOINTS ") ~= 0)
    n_datapoints = str2double(values{contains(params, "N_DATAPOINTS ")});
end
n_sites = max(mt_lengths);
n_dims = 2;

posFile = sprintf('%s/%s_filament_pos.file', file_dir, sim_name);
occuFile = sprintf('%s/%s_occupancy.file', file_dir, sim_name);
proteinFile = sprintf('%s/%s_protein_ids.file', file_dir, sim_name);

% filament position data; gives coordinates of each endpoint
filament_pos = zeros(n_mts, n_datapoints);
if isfile(posFile)
    mt_data_file = fopen(posFile);
    mt_data = fread(mt_data_file, 2*n_dims * n_mts * n_datapoints, '*double');
    fclose(mt_data_file);
    filament_pos = reshape(mt_data, n_dims, 2, n_mts, n_datapoints);
end
% occupancy data on each MT
occupancy = zeros(n_sites, n_mts, n_datapoints);
if isfile(occuFile)
    occupancy_file = fopen(occuFile);
    occu_raw_data = fread(occupancy_file, n_datapoints * n_mts * n_sites, '*int');
    occupancy = reshape(occu_raw_data, n_sites, n_mts, n_datapoints);
    fclose(occupancy_file);
end
% protein ID data
protein_ids = zeros(n_sites, n_mts, n_datapoints) - 1;
if isfile(proteinFile)
    protein_data_file = fopen(proteinFile);
    protein_raw_data = fread(protein_data_file, [n_mts * n_sites * n_datapoints], '*int');
    fclose(protein_data_file);
    protein_ids = reshape(protein_raw_data, n_sites, n_mts, n_datapoints);
end

site_ID = 0;
xlink_ID = 1;
motor_ID = 2;

motor_matrix = occupancy;
site_matrix = occupancy;

motor_matrix(motor_matrix ~= motor_ID) = 0;
motor_matrix(motor_matrix == motor_ID) = scale_factor;
site_matrix(site_matrix == xlink_ID) = -1;
site_matrix(site_matrix == motor_ID) = -1;
%site_matrix(site_matrix == site_ID) = 1;
site_matrix(site_matrix == -1) = 0;

if i_end == -1
    i_end = n_datapoints;
end

dwell_steps = dwell_time / time_per_datapoint;

pixels_x = ceil(n_sites*siteLength/pixelLength)+2*pixelPad;
if n_mts == 2
    pixels_x = ceil(2*n_sites*siteLength/pixelLength)+2*pixelPad;
    max_diff_x = 0;
    for i_data = i_start : dwell_steps : i_end - dwell_steps
        minus_one = filament_pos(:, 2, 1, i_data);
        plus_two = filament_pos(:, 1, 2, i_data);
        diff_x = plus_two(1) - minus_one(1);
        if diff_x > max_diff_x
           max_diff_x = diff_x; 
        end
    end
    n_sites_diff_max = ceil(max_diff_x / siteLength);
    pixel_diff_max = ceil(max_diff_x / pixelLength);
    pixels_x = pixels_x + (pixel_diff_max - 1);
end
pixels_y = ceil((i_end - i_start) / dwell_steps);
final_img = zeros(pixels_y, pixels_x, 3);
 % RGB image; 

% Run through movie frames and create each one
for i_data = i_start : dwell_steps :i_end - dwell_steps
    if n_mts == 2
        minus_one = filament_pos(:, 2, 1, i_data);
        plus_two = filament_pos(:, 1, 2, i_data);
        diff_x = plus_two(1) - minus_one(1);
        pixel_diff = ceil(diff_x / pixelLength);
        n_sites_diff = ceil(diff_x / siteLength);
        buffer = zeros(1, n_sites_diff);
        leftover = zeros(1, n_sites_diff_max - n_sites_diff);
        
        motors1 = sum(motor_matrix(:, 1, i_data:i_data + dwell_steps), 3)';
        motors2 = sum(motor_matrix(:, 2, i_data:i_data + dwell_steps), 3)';
        dataMatrix = [motors1 buffer motors2 leftover];
        
        sites1 = sum(site_matrix(:, 1, i_data:i_data + dwell_steps), 3)';
        sites2 = sum(site_matrix(:, 2, i_data:i_data + dwell_steps), 3)';
        lineMatrix = [sites1 buffer sites2 leftover];
    else
        dataMatrix = sum(motor_matrix(:, :, i_data:i_data + dwell_steps), 3)';
        lineMatrix = sum(site_matrix(:, :, i_data:i_data + dwell_steps), 3)';
    end
    
    % green channel - motors
    imageMotors = imageGaussianOverlap(dataMatrix,siteLength,pixelLength,pixelPad,...
        gaussSigma,gaussAmp,bkgLevel,noiseStd,doPlot);  
    imageMotors = imageMotors/intensityMax; %convert to grayscale
    % red channel - microtubule
    %lineMatrix = ones(size(dataMatrix));
   
    imageLine = imageGaussianOverlap(lineMatrix,siteLength,pixelLength,pixelPad,...
        gaussSigma,gaussAmp,bkgLevel,noiseStd,doPlot);
    imageLine = imageLine/intensityMax; %convert to grayscale
    % blue channel - noise
    imageBlue = ones(size(imageLine))*bkgLevel + randn(size(imageLine))*noiseStd;
    imageBlue = imageBlue/intensityMax;
    % merge into RGB image
    imageRGB = cat(3, imageLine, imageMotors, imageBlue);
    index = (i_data - i_start) / dwell_steps + 1;
    final_img(index, :, :) = mean(imageRGB(:, :, :), 1);
end
%} 
fig1 = figure;
set(fig1, 'Position', [100 100 500 700]);
axes('Units','Normalize','Position',[0 0 1 1]);
imagesc(final_img);
set(gca,'Xtick',[]); set(gca,'Ytick',[]);
set(gca, 'Box', 'off');
% Add a scale bar
len_x = scale_x * 1000 / pixelLength;
len_y = scale_t / dwell_time;
l = line([15 len_x+15],[5 5],'Color','w','LineWidth',4);
%l = line([9 len_x+9],[200 200],'Color','w','LineWidth',4);
set(l,'clipping','off')
l2 = line([10 10],[10 len_y+10],'Color','w','LineWidth',4);
%l2 = line([10 10],[350 len_y+350],'Color','w','LineWidth',4);
set(l2,'clipping','off')
