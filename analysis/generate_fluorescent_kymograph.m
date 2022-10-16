clear variables; 

file_dir = '/home/shane/projects/CyLaKS';
sim_name = 'shep_multiPF_1nM_20nM_0.0131_8';


dwell_time = 0.1;  % dwell time of theoretical camera
i_start = 1; %1150; %4150;
i_end = -1; %1650; %4650;
frac_visible = [1, 2]; %[2,3] % [numerator, denominator]; [1,1] for all visibile
xlink_intensity = (1/8); %1.5 % Controls how bright a single motor is (1 typically)
motor_intensity = (1/8);

subfilaments = true; 

% Scale bar lengths 
scale_x = 2; %2.5; %1; % microns
scale_t = 60; %30; %10; % seconds

% parameters for making simulated image (i.e., each frame)
siteLength = 8.2;
pixelLength = 150;
pixelPad = 5;
gaussSigma = 1.0;
doPlot = 0;
gaussAmp = 4000;
bkgLevel = 200;
noiseStd = 100;
intensityMax = gaussAmp/2 + bkgLevel;

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
if any(contains(params, "n_subfilaments") ~= 0)
    n_sub = sscanf(values{contains(params, "n_subfilaments ")}, '%g');
    if n_sub > n_mts
       n_mts = n_sub;
    end
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

n_dims = 2;

posFile = sprintf('%s/%s_filament_pos.file', file_dir, sim_name);
occuFile = sprintf('%s/%s_occupancy.file', file_dir, sim_name);
proteinFile = sprintf('%s/%s_protein_id.file', file_dir, sim_name);

% filament position data; gives coordinates of each endpoint
filament_pos = zeros(n_mts, n_datapoints);
if isfile(posFile)
    mt_data_file = fopen(posFile);
    mt_data = fread(mt_data_file, 2*n_dims * n_mts * n_datapoints, '*double');
    fclose(mt_data_file);
    filament_pos = reshape(mt_data, n_dims, 2, n_mts, n_datapoints);
end
% occupancy data on each MT
occupancy = zeros(max_sites, n_mts, n_datapoints);
if isfile(occuFile)
    occupancy_file = fopen(occuFile);
    occu_raw_data = fread(occupancy_file, n_datapoints * n_mts * max_sites, '*int');
    occupancy = reshape(occu_raw_data, max_sites, n_mts, n_datapoints);
    fclose(occupancy_file);
end
% protein ID data
protein_ids = zeros(max_sites, n_mts, n_datapoints) - 1;
if isfile(proteinFile)
    protein_data_file = fopen(proteinFile);
    protein_raw_data = fread(protein_data_file, [n_mts * max_sites * n_datapoints], '*int');
    fclose(protein_data_file);
    protein_ids = reshape(protein_raw_data, max_sites, n_mts, n_datapoints);
end

site_ID = 0;
xlink_ID = 1;
motor_ID = 2;

motor_matrix = zeros(mt_lengths(1), n_mts, n_datapoints); %protein_ids;
xlink_matrix = zeros(mt_lengths(1), n_mts, n_datapoints);
for i_data = 1 : n_datapoints
    for i_mt = 1 : n_mts
        for i_site = 1 : mt_lengths(i_mt)
            id = protein_ids(i_site, i_mt, i_data);
            sid = occupancy(i_site, i_mt, i_data);
            if sid == 1 %~= -1
                    xlink_matrix(i_site, i_mt, i_data) = xlink_intensity;
            elseif sid == 2
                %if mod(frac_visible(1)*id, frac_visible(2)) == 0
                    motor_matrix(i_site, i_mt, i_data) = motor_intensity;
                %end
            end
        end
    end
end
% change from 'zeros' to 'ones' to make MTs fluorescent 
site_matrix = zeros(mt_lengths(1), n_mts, n_datapoints); %occupancy;

if i_end == -1
    i_end = n_datapoints;
end

dwell_steps = int32(dwell_time / time_per_datapoint);

pixels_x = ceil(max_sites*siteLength/pixelLength)+2*pixelPad;
if n_mts == 2 && subfilaments == false
    pixels_x = ceil(2*max_sites*siteLength/pixelLength)+2*pixelPad;
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
% Sometimes you gotta add or subtract 1 here -- must be some kind of
% rounding error 
if n_mts == 2 && subfilaments == false
   pixels_x = pixels_x + 1; 
end
pixels_y = ceil((i_end - i_start) / dwell_steps);
final_img = zeros(pixels_y, pixels_x, 3);
 % RGB image; 
% Run through movie frames and create each one
for i_data = i_start : dwell_steps : i_end - dwell_steps
    if n_mts == 2 && subfilaments == false
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
        dataMatrixMotors = sum(motor_matrix(:, :, i_data:i_data + dwell_steps), [2 3])';
        dataMatrixXlinks = sum(xlink_matrix(:, :, i_data:i_data + dwell_steps), [2 3])';
        lineMatrix = sum(site_matrix(:, :, i_data:i_data + dwell_steps), [2 3])';
    end
    
    % green channel - motors
    imageMotors = imageGaussianOverlap(dataMatrixMotors,siteLength,pixelLength,pixelPad,...
        gaussSigma,gaussAmp,bkgLevel,noiseStd,doPlot);  
    imageMotors = imageMotors + ones(size(imageMotors))*bkgLevel + randn(size(imageMotors))*noiseStd;
    imageMotors = imageMotors/intensityMax; %convert to grayscale
    % red channel - microtubule
    %lineMatrix = ones(size(dataMatrix));
   
    imageLine = imageGaussianOverlap(lineMatrix,siteLength,pixelLength,pixelPad,...
        gaussSigma,gaussAmp,bkgLevel,noiseStd,doPlot);
    imageLine = imageLine/intensityMax; %convert to grayscale
    % blue channel - noise (now xlinks)
    
    imageBlue = ones(size(imageLine))*bkgLevel + randn(size(imageLine))*noiseStd;
    imageBlue = imageBlue/intensityMax;
    
    imageXlinks = imageGaussianOverlap(dataMatrixXlinks,siteLength,pixelLength,pixelPad,...
        gaussSigma,gaussAmp,bkgLevel,noiseStd,doPlot);  
    imageXlinks = imageXlinks + ones(size(imageXlinks))*bkgLevel + randn(size(imageXlinks))*noiseStd;
    imageXlinks = imageXlinks/intensityMax; %convert to grayscale
    % merge into RGB image
    imageRGB = cat(3, imageLine + imageMotors, imageXlinks, imageMotors);
    index = (i_data - i_start) / dwell_steps + 1;
    final_img(index, :, :) = mean(imageRGB(:, :, :), 1);
end
%} 
fig1 = figure;
set(fig1, 'Position', [100 100 350 350]);
%set(fig1, 'Position', [100 100 540 700]);
axes('Units','Normalize','Position',[0 0 1 1]);
img1 = imagesc(final_img);
set(gca,'Xtick',[]); set(gca,'Ytick',[]);
set(gca, 'Box', 'off');
% Add a scale bar
len_x = scale_x * 1000 / pixelLength;
len_y = scale_t / dwell_time;
%{
l = line([103.5 103.5-len_x],[2625 2625],'Color','w','LineWidth',4); %endtags
l2 = line([103 103],[2600 2600-len_y],'Color','w','LineWidth',4); %endtags
%}
%{
l = line([57 57-len_x],[980 980],'Color','w','LineWidth',4); %tubulin
l2 = line([56.75 56.75],[970 970-len_y],'Color','w','LineWidth',4); %tubulin
%}
x1 = (95/100)*pixels_x;
x2 = (94.25/100)*pixels_x;
y1 = (95/100)*pixels_y;
y2 = (93/100)*pixels_y;
l = line([x1 x1-len_x],[y1 y1],'Color','w','LineWidth',4); %tubulin
l2 = line([x2 x2],[y2 y2-len_y],'Color','w','LineWidth',4); %tubulin

set(l,'clipping','off')
set(l2,'clipping','off')
