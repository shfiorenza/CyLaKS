clear all;
% Often-changed variables
n_sites = 1000;
n_seeds = 24; 
initial_shift = [900,800,600,400,200,0];
n_shifts = length(initial_shift);
initial_overlap = (1000 - initial_shift)*0.008;
% Pseudo-constant variables
n_mts = 2;
n_steps = 100000000;
delta_t = 0.00001;
n_datapoints = 10000;
starting_point = 0;
% Calculate real time that passes per iteration to get an accurate velocity
time_per_datapoint = delta_t * (n_steps / n_datapoints);
% Calculate parameters for plotting / etc;
length = n_sites * 0.008;
end_time = n_steps * delta_t;
start_time = starting_point * delta_t;
fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
fileStructure = '%s_mt_coord.file';

max_velocities = zeros([n_seeds n_shifts]);
avg_max_velocities = zeros([1 n_shifts]);
err_max_velocities = zeros([1 n_shifts]);

max_velocities_b = zeros([n_seeds n_shifts]);
avg_max_velocities_b = zeros([1 n_shifts]);
err_max_velocities_b = zeros([1 n_shifts]);

final_lengths = zeros([n_seeds n_shifts]);
avg_final_lengths = zeros([1 n_shifts]);
err_final_lengths = zeros([1 n_shifts]);

for i_seed=1:n_seeds
    for i_shift=1:n_shifts
        simName = sprintf('scan_multi_slide/shift_%i_%i', ...
            i_seed-1, initial_shift(i_shift));
        
        simName_b = sprintf('scan_multi_stall_b/mt_coords/shift_stall_%i_%i', ...
            i_seed-1, initial_shift(i_shift));
        
        fileName = sprintf(fileDirectory, sprintf(fileStructure, simName));
        data_file = fopen(fileName);
        raw_data = fread(data_file, [n_mts, n_datapoints], 'double');
        fclose(data_file);
        fileName_b = sprintf(fileDirectory, sprintf(fileStructure, simName_b));
        data_file_b = fopen(fileName_b);
        raw_data_b = fread(data_file_b, [n_mts, n_datapoints], 'double');
        fclose(data_file_b);
        % Run through raw coord data to get overlap length at every datapoint
        final_overlap_data = zeros(n_datapoints, 1);
        final_overlap_data_b = zeros(n_datapoints, 1);
        for i_data=starting_point+1:n_datapoints
            mt_coord_one = raw_data(1, i_data);
            mt_coord_two = raw_data(2, i_data);
            delta = (mt_coord_two - mt_coord_one) * 0.008;
            overlap_length = abs(length - delta);
            final_overlap_data(i_data) = overlap_length;
            
            mt_coord_one = raw_data_b(1, i_data);
            mt_coord_two = raw_data_b(2, i_data);
            delta = (mt_coord_two - mt_coord_one) * 0.008;
            overlap_length = abs(length - delta);
            final_overlap_data_b(i_data) = overlap_length;
        end
        % Use gradient function with above spacing to get slope of overlap length
        final_overlap_data = smooth(final_overlap_data, 10);
        final_overlap_data_b = smooth(final_overlap_data_b, 10);
        %smoothed_overlap_data = final_overlap_data;
        slope_data = gradient(final_overlap_data, time_per_datapoint);
        slope_data_b = gradient(final_overlap_data_b, time_per_datapoint);
        final_velocity_data = smooth(slope_data, 100);
        final_velocity_data_b = smooth(slope_data_b, 100);
        % get max velocity during interval and convert to nm/s
        max_velocities(i_seed, i_shift) = max(abs(final_velocity_data)) * 1000; 
        max_velocities_b(i_seed, i_shift) = max(abs(final_velocity_data_b)) * 1000;
        % get average overlap length during last 1/4 of sim; should be steady
        final_lengths(i_seed, i_shift) = mean(final_overlap_data(0.75*n_datapoints:n_datapoints));
    end
end

for i_shift=1:n_shifts
    avg_final_lengths(i_shift) = mean(final_lengths(:, i_shift));
    avg_max_velocities(i_shift) = mean(max_velocities(:, i_shift));
    avg_max_velocities_b(i_shift) = mean(max_velocities_b(:, i_shift));
    % get standard devation of velocities
    err_max_velocities(i_shift) = std(max_velocities(:, i_shift));
    err_max_velocities_b(i_shift) = std(max_velocities_b(:, i_shift));
end

fig1 = figure();
set(fig1, 'Position', [50, 50, 1.8*480, 1.8*300])
%{
% Plot overlap length on top
subplot(2, 1, 1)
plot(initial_overlap, avg_final_lengths,'LineWidth', 2);
title(sprintf('Overlap length over time (%g microns or %d sites in length)', ...
        length, n_sites));
ylabel('Overlap length (microns)');
xlabel('Time (s)');
axis tight
%xlim([start_time + 1 end_time - 1]);
%ylim([0 n_sites * 0.008]);
grid on
grid minor
% Plot sliding velocity on bottom
subplot(2, 1, 2)
%}
errorbar(initial_overlap, avg_max_velocities, err_max_velocities, 'o-.','LineWidth',2);
hold on
%errorbar(initial_overlap, avg_max_velocities_b, err_max_velocities_b, 'o-.','LineWidth',2);
%title('Scaling of sliding velocity for 8-micron microtubules');
x = xlabel('Initial overlap length (microns)');
set(x, 'Units', 'Normalized', 'Position', [0.5, -0.065, 0]);
y = ylabel('Maximum sliding velocity (nm/s)');
set(y, 'Units', 'Normalized', 'Position', [-0.09, 0.5, 0]);
legend('0.5 nM PRC1 + 6 nM KIF4A','location', 'northwest', 'FontSize', 14);
xticks([0,2,4,6,8,10]);
yticks([0,100,200,300]);
xlim([0 9]);
ylim([0 350]);
axes = gca;
axes.XAxis.FontSize = 14;
axes.YAxis.FontSize = 14;
set(gcf, 'color', 'w');