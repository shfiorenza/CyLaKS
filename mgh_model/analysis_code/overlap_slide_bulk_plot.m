clear all;
% Effective binding concentration for tethered motors/xlinks
c_eff_teth = 4500;
% Effective binding concentration for second head of xlinks
c_eff_bind = 4500;
% Parameter to be scanned over for multiple scaling plots
scan_param = [10.0]; % c_xlink
legendLabel = cellstr(num2str(scan_param' / 10,'%.1f nM PRC1'));
% Length of top, mobile MT in n_sites
top_lengths = [75,150,225,300,375,450,525,600];
seeds = [1,2,3,4,5,6,7,8,9,11];
base_name = "2019_12_04_slideScan/slide_scan_%0.1f_%i_%i";
fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
site_size = 0.008;
%Experimental velocity data for 1 nM PRC1
exp_lengths_10 = [1.2, 1.6, 2.2, 2.7, 3.2, 4.0];
exp_vel_10 = [8, 13, 10, 10, 8, 5];
exp_err_10 = [4, 7, 4, 6, 9, 8];
%Experimental velocity data for 0.2 nM PRC1
exp_lengths_2 = [0.91, 1.5, 1.9, 2.1, 2.8, 4.0];
exp_vel_2 = [26, 40, 38, 38, 45, 61];
exp_err_2 = [18, 13, 32, 18, 12, 13];

final_lengths = zeros(length(scan_param), length(top_lengths), length(seeds));
peak_velocities = zeros(length(scan_param), length(top_lengths), length(seeds));

for i_scan=1:length(scan_param)
    for i_length=1:length(top_lengths)
        for i_seed=1:length(seeds)
            sim_name = sprintf(base_name, scan_param(i_scan), top_lengths(i_length), seeds(i_seed));
            log_file = sprintf(fileDirectory, sprintf('%s.log', sim_name));
            % Open log file and parse it into param labels & their values
            log = textscan(fileread(log_file),'%s %s', 'Delimiter', '=');
            params = log{1,1};
            values = log{1,2};
            % Read in number of MTs
            n_mts = str2double(values{contains(params, "count")});
            % Read in length of each MTs
            [length_one, length_two] = values{contains(params, "length")};
            mt_lengths(1) = sscanf(length_one, '%i');
            mt_lengths(2) = sscanf(length_two, '%i');
            % Read in auxiliary params
            n_steps = str2double(values{contains(params, "n_steps")});
            delta_t = sscanf(values{contains(params, "delta_t")}, '%g');
            n_datapoints = str2double(values{contains(params, "n_datapoints")});
            % Calculate parameters for plotting / etc;
            time_per_datapoint = delta_t * n_steps / n_datapoints;
            % Open mt coordinate file
            mt_coords_file = fopen(sprintf(fileDirectory, sprintf('%s_mt_coord.file', sim_name)));
            mt_coord_data = fread(mt_coords_file, [n_mts, n_datapoints], 'double');
            fclose(mt_coords_file);
            
            overlap_length = zeros(n_datapoints, 1);
            for i_data=1:n_datapoints
                mt_coords = mt_coord_data(:, i_data)';
                %{
                mt_endpoints = mt_coords + mt_lengths;
                overlap_start = max(mt_coords);
                overlap_end = min(mt_endpoints);
                %}
                plus_end_one = mt_coords(1);
                plus_end_two = mt_coords(2) + mt_lengths(2);
                
                %overlap_length(i_data) = (overlap_end - overlap_start) * site_size;
                overlap_length(i_data) = abs(plus_end_two - plus_end_one) * site_size;
            end
            % Get average final overlap length
            avg_length = mean(overlap_length(9333:10000));
            final_lengths(i_scan, i_length, i_seed) = avg_length;
            % Get max velocity during sliding
            overlap_length = smoothdata(overlap_length, 'movmean', 10);
            overlap_vel = smoothdata(gradient(overlap_length, time_per_datapoint), 'movmean', 100);
            peak_vel = max(abs(overlap_vel(6667:7333)));
            % Convert velocity from um/s to nm/s
            peak_velocities(i_scan, i_length, i_seed) = peak_vel * 1000;
        end
    end
end

avg_final_lengths = zeros(length(scan_param), length(top_lengths));
err_final_lengths = zeros(length(scan_param), length(top_lengths));
avg_velocities = zeros(length(scan_param), length(top_lengths));
err_velocities = zeros(length(scan_param), length(top_lengths));
for i_scan=1:length(scan_param)
    for i_length=1:length(top_lengths)
        % Get average final length across all seeds
        mean_final_length = mean(final_lengths(i_scan, i_length, :));
        avg_final_lengths(i_scan, i_length) = mean_final_length;
        % Calculate standard deviation of this average
        variance = 0.0;
        for i_seed=1:length(seeds)
            diff_sq = (mean_final_length - final_lengths(i_scan, i_length, i_seed))^2;
            variance = variance + diff_sq/(length(seeds) - 1);
        end
        err_final_lengths(i_scan, i_length) = sqrt(variance);
        % Get average sliding velocity across all seeds
        mean_vel = mean(peak_velocities(i_scan, i_length, :));
        avg_velocities(i_scan, i_length) = mean_vel;
        % Calculate standard deviation of this average
        variance = 0.0;
        for i_seed=1:length(seeds)
            diff_sq = (mean_vel - peak_velocities(i_scan, i_length, i_seed))^2;
            variance = variance + diff_sq/(length(seeds) - 1);
        end
        err_velocities(i_scan, i_length) = sqrt(variance);
    end
end

fig1 = figure();
set(fig1, 'Position', [50, 50, 1080, 720]);
hold all
for i_scan=1:length(scan_param)
    errorbar(top_lengths*site_size, avg_velocities(i_scan, :), ...
        err_velocities(i_scan, :), 'o', 'LineWidth', 2);
end
%Plot experimental line of best fit - length
%plot([1.0 4.5], [0.125 0.54], 'r', 'LineWidth', 2)
%legend({'Simulation data', 'Experimental line of best fit'}, 'Location', 'northwest');

%Plot experimental data - velocities
errorbar(exp_lengths_2, exp_vel_2, exp_err_2, 's', 'LineWidth', 2);
legend({'Simulation data for 1 nM PRC1', 'Experimental data for 0.2 nM PRC1'}, 'Location', 'northwest');

%{
for i_seed=1:length(seeds)
    plot(top_lengths*site_size, final_data(:, :, i_seed), '-o', 'LineWidth', 2);
end
%}

title('1 nM PRC1 + 6 nM Kif4A', 'FontSize', 14);
ylabel('Peak velocity during sliding (nm/s)', 'FontName', 'Garuda', 'FontSize', 12);
xlabel('Initial overlap length (microns)', 'FontName', 'Garuda', 'FontSize', 12)
xlim([0 5.5]);
axes = gca;
axes.XAxis.FontSize = 14;
axes.YAxis.FontSize = 14;
if length(scan_param) > 1
    legend(legendLabel, 'Location', 'best');
end