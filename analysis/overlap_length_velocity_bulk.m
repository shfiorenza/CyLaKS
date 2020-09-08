clear variables;
% Often-changed variables
n_sites = [1000,500];
n_sites_max = max(n_sites);
eff_concs = [1500,3000,4500];
k_unteth_scale = [1,5,25,125,625,3125,15625,78125];
base_name = 'slide_%i_%i';
folder_name = 'old_output/scan_output/%s';
n_concs = length(eff_concs);
n_rates = length(k_unteth_scale);
% Pseudo-constant variables
n_mts = length(n_sites);
n_steps = 100000000;
delta_t = 0.00001;
n_datapoints = 10000;
starting_point = 0;
% Calculate real time that passes per iteration to get an accurate velocity
time_per_datapoint = delta_t * (n_steps / n_datapoints);
% Calculate parameters for plotting / etc;
end_time = n_steps * delta_t;
start_time = starting_point * delta_t;
% Directory where files are located
fileDirectory = sprintf('/home/shane/Desktop/%s', folder_name);
%fileDirectory = sprintf('/home/shane/Projects/overlap_analysis/mgh_model/%s', ...
%    folder_name);
fileStructure = '%s_mt_coord.file';

max_velocities = zeros([n_concs n_rates]);
for i_conc=1:n_concs
    for i_rate=1:n_rates
        simName = sprintf(base_name, eff_concs(i_conc), k_unteth_scale(i_rate));
        fileName = sprintf(fileDirectory, sprintf(fileStructure, simName));
        data_file = fopen(fileName);
        data_exists = true;
        try 
            raw_data = fread(data_file, [n_mts, n_datapoints], 'double');
        catch
            disp('wut')
            max_velocities(i_conc, i_rate) = -5;
            data_exists = false;
        end
        fclose(data_file);
        if data_exists
            final_overlap_data = zeros(n_datapoints, 1);
            % Run through raw coord data to get overlap length at every datapoint
            for i_data=starting_point+1:n_datapoints
                mt_coord_one = raw_data(1, i_data);
                mt_coord_two = raw_data(2, i_data);
                plus_end_one = mt_coord_one;
                plus_end_two = mt_coord_two + n_sites(2);
                end_dist = (plus_end_two - plus_end_one) * 0.008;
                final_overlap_data(i_data) = end_dist;
            end
            % Use gradient function with above spacing to get slope of overlap length
            final_overlap_data = smoothdata(final_overlap_data, 'movmean',300);
            slope_data = gradient(final_overlap_data, time_per_datapoint);
            final_velocity_data = smoothdata(slope_data, 'movmean', 300);
            % Get max velocity during interval and convert to nm/s
            max_velocities(i_conc, i_rate) = max(abs(final_velocity_data(4500:5500))) * 1000;
            
        end
    end
end
k_unteth = k_unteth_scale * 0.0005;
fig1 = figure();
set(fig1, 'Position', [50, 50, 1.8*480, 1.8*300])
for i_conc=1:n_concs
   semilogx(k_unteth, max_velocities(i_conc, :), '*', 'MarkerSize', 16, 'LineWidth', 2);
   hold all
end
title('Antiparallel overlap length of 4 microns', 'FontSize', 14);
xlabel({'Untethering rate (s^{-1})','Tethering kD held constant at 30 nM'}, 'FontSize', 14);
ylabel('Maximum sliding velocity (nm/s)', 'FontSize', 14);
ylim([85 120]);
xlim([1e-4 2e2]);
legendLabel = cellstr(num2str(eff_concs', '%i nM effective tether conc.'));
legend(legendLabel, 'location', 'northwest', 'FontSize', 13);
axes = gca;
axes.XAxis.FontSize = 14;
axes.YAxis.FontSize = 14;
yticks([90,100, 110, 120]);
%errorbar(initial_overlap, avg_max_velocities, err_max_velocities, 'o-.','LineWidth',2);
%errorbar(initial_overlap, avg_max_velocities_b, err_max_velocities_b, 'o-.','LineWidth',2);
%title('Scaling of sliding velocity for 8-micron microtubules');
%{
x = xlabel('Initial overlap length (microns)');
set(x, 'Units', 'Normalized', 'Position', [0.5, -0.065, 0]);
y = ylabel('Maximum sliding velocity (nm/s)');
set(y, 'Units', 'Normalized', 'Position', [-0.09, 0.5, 0]);
legend('0.5 nM PRC1 + 6 nM KIF4A','location', 'northwest', 'FontSize', 14);

xticks([0,1,2,3]); %,6,8,10]);
%yticks([0,10,20,30,40,50]);
%xlim([0 3]);
%ylim([0 55]);
axes = gca;
axes.XAxis.FontSize = 14;
axes.YAxis.FontSize = 14;
set(gcf, 'color', 'w');
%}