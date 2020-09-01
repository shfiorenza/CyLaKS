
clear all;
% Often-changed variables
n_sites = 1000;
k_D = [300,100,30];%,15,3]; 
suffix = ["e","i","d"];%,"g","c"];
n_runs = length(suffix);
simName = 'slide_mucho_teth_%s';
% Pseudo-constant variables
n_mts = 2;
n_steps = 100000000;
delta_t = 0.00001;
n_datapoints = 10000;
% Calculate parameters for plotting / etc;
length = n_sites * 0.008;
end_time = n_steps * delta_t;

fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
fileStructure = '%s_mt_coord.file';

% Run through raw coord data to get overlap length at every datapoint
final_overlap_data = zeros([n_runs n_datapoints]);
for i_run=1:n_runs
    for i_data = 1:n_datapoints
        fileName = sprintf(fileDirectory, sprintf(fileStructure, ...
            sprintf(simName, suffix(i_run))));
        data_file = fopen(fileName);
        raw_data = fread(data_file, [n_mts, n_datapoints], 'double');
        fclose(data_file);
        
        mt_coord_one = raw_data(1, i_data);
        mt_coord_two = raw_data(2, i_data);
        delta = (mt_coord_two - mt_coord_one) * 0.008;
        if(delta > 0)
            overlap_length = length - delta;
        else
            overlap_length = length + delta;
        end
        final_overlap_data(i_run, i_data) = overlap_length;
    end
    final_overlap_data(i_run, :) = smooth(final_overlap_data(i_run, :), 50);
end

fig1 = figure();
set(fig1, 'Position', [50, 50, 3*480, 2*300])
plot(linspace(-500, end_time-500, n_datapoints), final_overlap_data, ...
    'LineWidth', 3);
yline(0, '--', 'LineWidth', 1);
%title(sprintf('Sliding of two 8-micron MTs for different tether dissociation constants', ...
%   length, n_sites),'FontSize', 17);
ylabel('Overlap length (microns)');
xlabel('Time (s)');
label = {'k_D = 300 nM', 'k_D = 100 nM','k_D = 30 nM'};%, 'k_D = 15 nM', 'k_D = 3 nM'};
legend(label, 'location', 'northeast', 'FontSize', 14);
xticks([0,100,200,300]);
yticks([-8,-4,0,4,8]);
xlim([-10 (end_time-150)-500]);
ylim([-1 9]);
axes = gca;
axes.XAxis.FontSize = 14;
axes.YAxis.FontSize = 14;
background = gcf;
set(background, 'color', 'w');