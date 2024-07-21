clear variables; 

%sim_name = 'out_endtags1/shep_1nM_100nM_8_250_0.2kT_1x_0';
sim_name = 'output16/shep_0.1nM_10nM_8_0.2kT_1x_0';
sim_name = 'motility_75_40_50nM_1nM_8x';
sim_name = 'output18/shep_1nM_100nM_8_1000_0.8kT_1x_0';

dwell_time = 0.1;  % dwell time of theoretical camera
i_start = 1;
i_end = -1; 
frac_visible_xlink = [1, 1]; % [numerator, denominator]; [1,1] for all visibile
frac_visible_motor = [1, 1]; % [numerator, denominator]; [1,1] for all visibile

tubulin_intensity = 0.0; % 0.01;
xlink_intensity = 0.0025; % Controls how bright a single xlink is 
motor_intensity = 0.01; %0.006;
subfilaments = true; 

% Scale bar lengths 
scale_x = 1; %2.5; %1; % microns
scale_t = 60; %30; %10; % seconds
% parameters for making simulated image (i.e., each frame)
siteLength = 8.2;
pixelLength = 150;
pixelPad = 2;
gaussSigma = 1.0;
doPlot = 0;
gaussAmp = 40000;
bkgLevel = 200;
noiseStd = 100;
intensityMax = gaussAmp/2 + bkgLevel;

% Load parameter structure
file_dir = '..';  % Default; only change if you move CyLaKS output files
params = load_parameters(sprintf('%s/%s', file_dir, sim_name));

% Open data files 
filament_filename = sprintf('%s/%s_filament_pos.file', file_dir, sim_name);
filament_pos = zeros(params.n_dims, 2, params.n_mts, params.n_datapoints);
filament_pos = load_data(filament_pos, filament_filename, '*double'); 

protein_filename = sprintf('%s/%s_protein_id.file', file_dir, sim_name);
protein_ids = zeros(params.max_sites, params.n_mts, params.n_datapoints) - 1;
protein_ids = load_data(protein_ids, protein_filename, '*int');

occupancy_filename = sprintf('%s/%s_occupancy.file', file_dir, sim_name);
occupancy = zeros(params.max_sites, params.n_mts, params.n_datapoints) - 1;
occupancy = load_data(occupancy, occupancy_filename, '*int');

site_ID = 0;
xlink_ID = 1;
motor_ID = 2;
if i_end == -1
    i_end = params.n_datapoints;
end
if subfilaments % rescale so we can compare across different n_subfilaments
    xlink_intensity = xlink_intensity/sqrt(params.n_mts);
    motor_intensity = motor_intensity/sqrt(params.n_mts);
end
dwell_steps = int32(dwell_time / params.time_per_datapoint);

% change from 'zeros' to 'ones' to make MTs fluorescent 
site_matrix = zeros(params.max_sites, params.n_mts, params.n_datapoints);
% Get motor and xlink matrices from occupancy data
motor_matrix = zeros(params.max_sites, params.n_mts, params.n_datapoints); % from protein_ids;
xlink_matrix = zeros(params.max_sites, params.n_mts, params.n_datapoints); % from protein_ids;
n_xlinks_avg = 0;
for i_data = 1 : params.n_datapoints
    for i_mt = 1 : params.n_mts
        for i_site = 1 : params.mt_lengths(i_mt)
            id = protein_ids(i_site, i_mt, i_data); % unique individual ID
            sid = occupancy(i_site, i_mt, i_data);  % species label ID 
            if sid == xlink_ID
                if mod(frac_visible_xlink(1)*id, frac_visible_xlink(2)) == 0
                    xlink_matrix(i_site, i_mt, i_data) = xlink_intensity;
                    n_xlinks_avg = n_xlinks_avg + 1/(params.n_datapoints);
                end
            elseif sid == motor_ID
                if mod(frac_visible_motor(1)*id, frac_visible_motor(2)) == 0
                    motor_matrix(i_site, i_mt, i_data) = motor_intensity;
                end
            end
        end
    end
end

disp(n_xlinks_avg)

% Calculate final image dimensions
pixels_y = ceil((i_end - i_start) / dwell_steps);
pixels_x = ceil(params.max_sites*siteLength/pixelLength)+2*pixelPad;
if params.n_mts == 2 && subfilaments == false
    max_span = 0;
    for i_data = i_start : dwell_steps : i_end - dwell_steps
        plus_one = filament_pos(1, 2, 1, i_data); % actually minus
        plus_two = filament_pos(1, 2, 2, i_data); % actually minus
        span = max(plus_one, plus_two) - min(plus_one, plus_two);
        if span > max_span
           max_span = span;
        end
    end
    span_sites_max = int32(max_span / siteLength);
    pixels_x = ceil(max_span/pixelLength) + 2*pixelPad;
    %{
    pixels_x = ceil(2*params.max_sites*siteLength/pixelLength)+2*pixelPad;
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
    %}
end

final_img_motors = zeros(pixels_y, pixels_x, 3); % RGB image; 
final_img_xlinks = zeros(pixels_y, pixels_x, 3); % RGB image; 
final_img_combined = zeros(pixels_y, pixels_x, 3); % RGB image; 
% Run through data and create each line of kymograph step-by-step
for i_data = i_start : dwell_steps : i_end - dwell_steps
    if params.n_mts == 2 && subfilaments == false
        % for sliding 
        plus_one = filament_pos(1, 2, 1, i_data); % actually minus
        plus_two = filament_pos(1, 2, 2, i_data); % actually minus
        span = max(plus_one, plus_two) - min(plus_one, plus_two);
        span_sites = int32(span / siteLength);
        i_offset = span_sites - params.mt_lengths(2);
        if i_offset < 1
            i_offset = 1;
        end
        
        xlinks1 = sum(xlink_matrix(:, 2, i_data:i_data + dwell_steps), 3)';
        xlinks2 = sum(xlink_matrix(:, 1, i_data:i_data + dwell_steps), 3)';
        dataMatrixXlinks = zeros(1, span_sites);
        dataMatrixXlinks(1:length(xlinks1)) = xlinks1; 
        dataMatrixXlinks(i_offset:i_offset+length(xlinks2)-1) = dataMatrixXlinks(i_offset:i_offset+length(xlinks2)-1) + xlinks2;
        dataMatrixXlinks(i_offset:length(xlinks2)-1) = dataMatrixXlinks(i_offset:length(xlinks2)-1) / 2;
        leftover = zeros(1, (span_sites_max - length(dataMatrixXlinks))/2);
        dataMatrixXlinks = [leftover dataMatrixXlinks leftover];
       % dataMatrixXlinks = [dataMatrixXlinks leftover leftover];
        
        lineMatrix = zeros(1, span_sites);
        lineMatrix(1:length(xlinks1)) = tubulin_intensity * ones(1, length(xlinks1)); 
        lineMatrix(i_offset:i_offset+length(xlinks2)-1) = lineMatrix(i_offset:i_offset+length(xlinks2)-1) + tubulin_intensity * ones(1, length(xlinks2)); 
        lineMatrix = [leftover lineMatrix leftover];
       % lineMatrix = [lineMatrix leftover leftover];
        
        dataMatrixMotors = zeros(1,length(dataMatrixXlinks));
        % for ablation
        %{
        minus_one = filament_pos(:, 2, 1, i_data);
        plus_two = filament_pos(:, 1, 2, i_data);
        diff_x = plus_two(1) - minus_one(1);
        pixel_diff = ceil(diff_x / pixelLength);
        n_sites_diff = ceil(diff_x / siteLength);
        buffer = zeros(1, n_sites_diff);
        leftover = zeros(1, n_sites_diff_max - n_sites_diff);
        
        motors1 = sum(motor_matrix(:, 1, i_data:i_data + dwell_steps), 3)';
        motors2 = sum(motor_matrix(:, 2, i_data:i_data + dwell_steps), 3)';
        dataMatrixMotors = [motors1 buffer motors2 leftover];
        
        xlinks1 = sum(xlink_matrix(:, 1, i_data:i_data + dwell_steps), 3)';
        xlinks2 = sum(xlink_matrix(:, 2, i_data:i_data + dwell_steps), 3)';
        dataMatrixXlinks = [xlinks1 xlinks2 buffer leftover];
        
        sites1 = sum(site_matrix(:, 1, i_data:i_data + dwell_steps), 3)';
        sites2 = sum(site_matrix(:, 2, i_data:i_data + dwell_steps), 3)';
        lineMatrix = [sites1 sites2 buffer leftover];
        %}
    else
        dataMatrixMotors = sum(motor_matrix(:, :, i_data:i_data + dwell_steps), [2 3])';
        dataMatrixXlinks = sum(xlink_matrix(:, :, i_data:i_data + dwell_steps), [2 3])';
        lineMatrix = sum(site_matrix(:, :, i_data:i_data + dwell_steps), [2 3])';
    end
    
    % Motors - purple (red + blue) channel
    imageMotors = imageGaussianOverlapSlice(dataMatrixMotors,siteLength,pixelLength,pixelPad,...
        gaussSigma,gaussAmp,bkgLevel,noiseStd,doPlot);
    imageMotors = imageMotors + ones(size(imageMotors))*bkgLevel + randn(size(imageMotors))*noiseStd;
    imageMotors = imageMotors/intensityMax;
    
    % Microtubules - red channel
    %lineMatrix = ones(size(dataMatrix));
    imageLine = imageGaussianOverlapSlice(lineMatrix,siteLength,pixelLength,pixelPad,...
        gaussSigma,gaussAmp,bkgLevel,noiseStd,doPlot);
    imageLine = imageLine/intensityMax; %convert to grayscale
    %imageLine = zeros(size(imageMotors));
    
    % Crosslinkers - green channel 
    imageXlinks = imageGaussianOverlapSlice(dataMatrixXlinks,siteLength,pixelLength,pixelPad,...
        gaussSigma,gaussAmp,bkgLevel,noiseStd,doPlot);  
    imageXlinks = imageXlinks + ones(size(imageXlinks))*bkgLevel + randn(size(imageXlinks))*noiseStd;
    imageXlinks = imageXlinks/intensityMax;

    % merge into RGB image
    imageNull = zeros(size(imageLine));

    imageRGB_motor = cat(3, imageMotors, imageNull, imageMotors);
    imageRGB_xlink = cat(3, imageNull, imageXlinks, imageNull);
    imageRGB_combined = cat(3, imageLine + imageMotors, imageXlinks, imageMotors);
    
    index = (i_data - i_start) / dwell_steps + 1;

    final_img_motors(index, :, :) = imageRGB_motor;
    final_img_xlinks(index, :, :) = imageRGB_xlink;
    final_img_combined(index, :, :) = imageRGB_combined;
end
%} 
fig_motor = figure;
set(fig_motor, 'Position', [50 50 300 600]);
axes('Units','Normalize','Position',[0 0 1 1]);
img1 = imagesc(final_img_motors, [min(final_img_motors, [], 'all') max(final_img_motors, [], 'all')]);
set(gca,'Xtick',[]); set(gca,'Ytick',[]);
set(gca, 'Box', 'off');
% Add scale bars
len_x = scale_x * 1000 / pixelLength;
len_y = scale_t / dwell_time;
x1 = (95/100)*pixels_x;
x2 = (94.25/100)*pixels_x;
y1 = (99/100)*pixels_y;
y2 = (97/100)*pixels_y;
l = line([x1 x1-len_x],[y1 y1],'Color','w','LineWidth',4); %tubulin
l2 = line([x2 x2],[y2 y2-len_y],'Color','w','LineWidth',4); %tubulin
set(l,'clipping','off')
set(l2,'clipping','off')

fig_xlink = figure;
set(fig_xlink, 'Position', [100 50 300 600]);
axes('Units','Normalize','Position',[0 0 1 1]);
img2 = imagesc(final_img_xlinks, [min(final_img_xlinks, [], 'all') max(final_img_xlinks, [], 'all')]);
set(gca,'Xtick',[]); set(gca,'Ytick',[]);
set(gca, 'Box', 'off');
% Add scale bars
len_x = scale_x * 1000 / pixelLength;
len_y = scale_t / dwell_time;
x1 = (95/100)*pixels_x;
x2 = (94.25/100)*pixels_x;
y1 = (99/100)*pixels_y;
y2 = (97/100)*pixels_y;
l = line([x1 x1-len_x],[y1 y1],'Color','w','LineWidth',4); %tubulin
l2 = line([x2 x2],[y2 y2-len_y],'Color','w','LineWidth',4); %tubulin
set(l,'clipping','off')
set(l2,'clipping','off')

fig_combined = figure;
set(fig_combined, 'Position', [150 50 300 600]);
axes('Units','Normalize','Position',[0 0 1 1]);
img3 = imagesc(final_img_combined, [min(final_img_combined, [], 'all') max(final_img_combined, [], 'all')]);
set(gca,'Xtick',[]); set(gca,'Ytick',[]);
set(gca, 'Box', 'off');
% Add scale bars
len_x = scale_x * 1000 / pixelLength;
len_y = scale_t / dwell_time;
x1 = (95/100)*pixels_x;
x2 = (94.25/100)*pixels_x;
y1 = (99/100)*pixels_y;
y2 = (97/100)*pixels_y;
l = line([x1 x1-len_x],[y1 y1],'Color','w','LineWidth',4); %tubulin
l2 = line([x2 x2],[y2 y2-len_y],'Color','w','LineWidth',4); %tubulin
set(l,'clipping','off')
set(l2,'clipping','off')

%saveas(fig_motor, 'kymo_motors.png', 'png');
%saveas(fig_xlink, 'kymo_xlinks.png', 'png');
%saveas(fig_combined, 'kymo_combo.png', 'png');
