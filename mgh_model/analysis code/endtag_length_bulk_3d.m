
clear variables;
off_ratio = 10;
%k_hydrolyzes = [60,75,90,105,120];
k_hydrolyzes = [60,90,120];
%jam_ratios = [100,200,300,400,500,600,700,800];
jam_ratios = [100,300,500,700];
% Pseudo-constant variables
mt_lengths = [2, 4, 6, 8, 10, 14]; % in microns
experiment_endtags = [1.1, 1.2, 1.35, 1.55, 1.8, 2.2];
exp_endtags = [1.15,1.25,1.45,1.6,1.8,2.15]; 
exp_errors = [0.1,0.1,0.2,0.3,0.25,0.4];
motor_speciesID = 2;
xlink_speciesID = 1;
n_datapoints = 10000;
starting_point = 5000;
active_datapoints = n_datapoints - starting_point;
fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
fileStruct = '%s_occupancy.file';

% Set color order so lines have consistent colors for up to 8 k_off_pseudo values
ColorArray = {'orange', 'cyan', 'blue', 'green', 'yellow', ... 
    'dark grey','magenta', 'red'};
RGBArray = [0.91 0.41 0.17; 0.301 0.745 0.933; 0 0 1; 0 0.5 0; ...
    0.929 0.6940 0.125; 0.3 0.3 0.3; 1 0 1; 1 0 0];

n_lengths = length(mt_lengths);
n_hydros = length(k_hydrolyzes);
n_ratios = length(jam_ratios); 
final_data = zeros([n_ratios n_hydros n_lengths]);

for i_ratio=1:1:n_ratios
    jam_ratio = jam_ratios(i_ratio); 
    for i_hydro=1:1:n_hydros
        k_hydro = k_hydrolyzes(i_hydro);
        for i_length=1:1:n_lengths
            n_sites = mt_lengths(i_length) * 125;
            if(off_ratio == 1)
                simName = sprintf('scan_output/k_hydrolyze_%i/jam_ratio_%i/endtag_%i_%i_%i', ...
                    k_hydro, jam_ratio, k_hydro, jam_ratio, n_sites);
            else
                simName = sprintf('scan_output_b/off_ratio_%i/k_hydrolyze_%i/jam_ratio_%i/endtag_%i_%i_%i_%i', ...
                    off_ratio, k_hydro, jam_ratio, off_ratio, k_hydro, jam_ratio, n_sites);
            %    simName = sprintf('endtag_%i_%i_%i_%i_long', off_ratio, k_hydro, jam_ratio, n_sites);
            end
            fileName = sprintf(fileDirectory, sprintf(fileStruct, simName));
            data_file = fopen(fileName);
            motor_raw_data = fread(data_file, [n_sites, n_datapoints], '*int');
            xlink_raw_data = motor_raw_data;
            fclose(data_file);
            
            motor_avg_occupancy = zeros([n_sites 1]);
            xlink_avg_occupancy = zeros([n_sites 1]);
            
            motor_raw_data(motor_raw_data ~= motor_speciesID) = 0;
            motor_raw_data(motor_raw_data == motor_speciesID) = 1;
            
            xlink_raw_data(xlink_raw_data ~= xlink_speciesID) = 0;
            xlink_raw_data(xlink_raw_data == xlink_speciesID) = 1;
            
            % Read in and average occupancy data over all datapoints
            for i=starting_point:1:n_datapoints
                motor_avg_occupancy(:,1) = motor_avg_occupancy(:,1) + double(motor_raw_data(:,i))./active_datapoints;
                xlink_avg_occupancy(:,1) = xlink_avg_occupancy(:,1) + double(xlink_raw_data(:,i))./active_datapoints;
            end
            
            motor_smoothed_avg = smoothdata(motor_avg_occupancy);
            xlink_smoothed_avg = smoothdata(xlink_avg_occupancy);
            net_smoothed_avg = motor_smoothed_avg + xlink_smoothed_avg;
            occupancy_slope = abs(gradient(net_smoothed_avg, 0.08));
            max_slope = max(occupancy_slope);
            
            endtag_site = 0;
            past_max = false;
            for i_site=1:n_sites
                if(occupancy_slope(i_site) >= max_slope)
                    past_max = true;
                end
                if(occupancy_slope(i_site) < max_slope/2 && past_max)
                    endtag_site = i_site;
                    break;
                end
            end
            endtag_length = endtag_site*0.008;
            final_data(i_ratio, i_hydro, i_length) = endtag_length;
        end
    end
end
%}



% x-axis always corresponds to MT length values
x_data = mt_lengths;
% y-axis corresponds to different jam ratios
for i_ratio=1:1:n_ratios
    % Make sure each jam_ratio "slice" has a constant y-value
    y_data(1:n_lengths) = jam_ratios(i_ratio);
    % z-axis will be avg occupancy 
    for i_hydro=1:1:n_hydros
        z_data = permute(final_data(i_ratio,i_hydro,:), [3 2 1]);
        plot3(x_data, y_data, z_data, '-*', 'Color', RGBArray(i_hydro, :), ...
            'LineWidth', 2);
        hold on
    end
end

% Style stuff
view(-77.5, 5)
set(gcf, 'Position', [50 50 1080 800]);
legendLabel = cellstr(num2str(k_hydrolyzes', 'k hydrolyze = %-g'));
legend(legendLabel, 'location', 'northeast', 'FontSize', 14, 'AutoUpdate', 'Off');
title(sprintf('Endtag scaling for k off ratio = %g', off_ratio));%, ...
xlabel('MT length (microns)');
xlabel_pos = get(get(gca, 'XLabel'), 'Position');
set(get(gca, 'XLabel'), 'Position', xlabel_pos + [0 0 -0.2]);
xticks(mt_lengths);
ylabel('Hydrolysis jam ratio');
ylabel_pos = get(get(gca, 'YLabel'), 'Position');
set(get(gca, 'YLabel'), 'Position', ylabel_pos + [0 0 -0.2]);
yticks(jam_ratios);
zlabel('Endtag length (microns)');
zlim([0 5]);

for i_ratio=1:1:n_ratios
   y_data(1:n_lengths) = jam_ratios(i_ratio);
   for i_hydro=1:1:n_hydros
        plot3(x_data, y_data, exp_endtags, '--r', 'LineWidth', 2);
    end
end
