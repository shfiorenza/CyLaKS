clear all;
% Often-changed variables
n_sites = [1000, 500];
%simName = sprintf('slide_scans/quart_proc_half_speed/slide_0_%i', n_sites(2));
simName = 'k_on_scan/slide_kOn_1000_10';
% Pseudo-constant variables
n_mts = 2;
n_steps = 100000000;
delta_t = 0.00001;
n_datapoints = 10000;
starting_point = 0;
% Calculate parameters for plotting / etc;
length = n_sites(2) * 0.008;
end_time = n_steps * delta_t;
start_time = starting_point * delta_t;

%fileDirectory = '/home/shane/Desktop/slide_scan/%s';
fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
fileStructure = '%s_mt_coord.file';
fileName = sprintf(fileDirectory, sprintf(fileStructure, simName));

data_file = fopen(fileName);
raw_data = fread(data_file, [n_mts, n_datapoints], 'double');
fclose(data_file);

% Run through raw coord data to get overlap length at every datapoint
final_overlap_data = zeros(n_datapoints, 1);
end_dist_data = zeros(n_datapoints, 1);
for i_data = starting_point + 1 : 1 : n_datapoints
    mt_coord_one = raw_data(1, i_data);
    mt_coord_two = raw_data(2, i_data);
    
    if(mt_coord_two + n_sites(2) <= mt_coord_one + n_sites(1) ...
    && mt_coord_two >= mt_coord_one)
        overlap_length = n_sites(2) * 0.008;
    else
        if(mt_coord_two > mt_coord_one)
           delta = (mt_coord_two + n_sites(2)) - (mt_coord_one + n_sites(1)); 
        else
           delta = mt_coord_one - mt_coord_two;
        end
        overlap_length = (n_sites(2) - delta) * 0.008;
    end
    final_overlap_data(i_data) = overlap_length; 
    
    plus_end_one = mt_coord_one;
    plus_end_two = mt_coord_two + n_sites(2);
    end_dist = (plus_end_two - plus_end_one) * 0.008;
    end_dist_data(i_data) = end_dist;
end

% Calculate real time that passes per iteration to get an accurate velocity
time_per_datapoint = delta_t * (n_steps / n_datapoints);
% Use gradient function with above spacing to get slope of overlap length
overlap_data = smoothdata(end_dist_data, 'movmean', 125);
slope_data = smoothdata(gradient(overlap_data, time_per_datapoint));

fig1 = figure();
set(fig1, 'Position', [50, 50, 2.5*480, 2.5*300])
hold on

% Plot overlap length on top
subplot(2, 1, 1)
plot(linspace(start_time, end_time, n_datapoints), overlap_data, ...
        'LineWidth', 2);
%title(sprintf('Overlap length over time (%g microns or %d sites in length)', ...
 %length, n_sites(2)));
 title('c eff teth = 1000 nM; k on xlink = 0.1 (nM * s)^{-1}');
ylabel('Overlap length (microns)');
xlabel('Time (s)');
axis tight
xlim([start_time end_time]);
ylim([0 n_sites(1) * 0.008]);
grid on
grid minor

% Plot sliding velocity on bottom
subplot(2, 1, 2)
plot(linspace(start_time, end_time, n_datapoints), slope_data * 1000, ...
        'LineWidth', 2);
title('Sliding velocity over time');
ylabel('Sliding velocity (nm/s)');
xlabel('Time (s)');
xlim([start_time end_time]);
%ylim([-250 250]);
%axis tight
grid on
grid minor

set(gcf, 'color', 'w');