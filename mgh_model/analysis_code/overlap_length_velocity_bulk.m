clear all;
% Often-changed variables
n_seeds = 1; %24;
bot_length = 1000; 
base_name = 'slide_%i_%i';
folder_name = 'slide_scans/quart_proc_half_speed/%s';
top_length = [100,150,200,250,300,350];
%top_length = [100,200,400,600,800,1000];
initial_overlap = top_length * 0.008;
n_lengths = length(top_length);
% Pseudo-constant variables
n_mts = 2;
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
fileDirectory = sprintf('/home/shane/Projects/overlap_analysis/mgh_model/%s', ...
    folder_name);
fileStructure = '%s_mt_coord.file';


max_velocities = zeros([n_seeds n_lengths]);
avg_max_velocities = zeros([1 n_lengths]);
err_max_velocities = zeros([1 n_lengths]);
%{
max_velocities_b = zeros([n_seeds n_shifts]);
avg_max_velocities_b = zeros([1 n_shifts]);
err_max_velocities_b = zeros([1 n_shifts]);
%}
final_lengths = zeros([n_seeds n_lengths]);
avg_final_lengths = zeros([1 n_lengths]);
err_final_lengths = zeros([1 n_lengths]);

for i_seed=1:n_seeds
    for i_length=1:n_lengths
        simName = sprintf(base_name, i_seed-1, top_length(i_length));
        fileName = sprintf(fileDirectory, sprintf(fileStructure, simName));
        data_file = fopen(fileName);
        raw_data = fread(data_file, [n_mts, n_datapoints], 'double');
        fclose(data_file);        
        final_overlap_data = zeros(n_datapoints, 1);     
        
        %{
        simName_b = sprintf('shift_%i_%i', i_seed-1, initial_shift(i_shift));
        fileName_b = sprintf(fileDirectory, sprintf(fileStructure, simName_b));
        data_file_b = fopen(fileName_b);
        raw_data_b = fread(data_file_b, [n_mts, n_datapoints], 'double');
        fclose(data_file_b); 
        final_overlap_data_b = zeros(n_datapoints, 1);      
        %}
        n_sites = [bot_length, top_length(i_length)];
        % Run through raw coord data to get overlap length at every datapoint
        for i_data=starting_point+1:n_datapoints
            mt_coords = raw_data(:, i_data);     

            if(mt_coords(2) + n_sites(2) <= mt_coords(1) + n_sites(1) ...
                    && mt_coords(2) >= mt_coords(1))
                overlap_length = n_sites(2) * 0.008;
            else
                if(mt_coords(2) > mt_coords(1))
                    delta = (mt_coords(2) + n_sites(2)) - (mt_coords(1) + n_sites(1));
                else
                    delta = mt_coords(1) - mt_coords(2);
                end
                overlap_length = (n_sites(2) - delta) * 0.008;
            end
            final_overlap_data(i_data) = overlap_length;
           
            %{
            mt_coord_one = raw_data_b(1, i_data);
            mt_coord_two = raw_data_b(2, i_data);
            delta = (mt_coord_two - mt_coord_one) * 0.008;
            overlap_length = abs(length - delta);
            final_overlap_data_b(i_data) = overlap_length;
            %}
        end     
        % Use gradient function with above spacing to get slope of overlap length
        final_overlap_data = smooth(final_overlap_data, 100);
        slope_data = gradient(final_overlap_data, time_per_datapoint);
        final_velocity_data = smooth(slope_data, 100);    
        % Get max velocity during interval and convert to nm/s      
        max_velocities(i_seed, i_length) = max(abs(final_velocity_data(1000:9000))) * 1000;
        % Get average overlap length during last 1/4 of sim; should be steady       
        final_lengths(i_seed, i_length) = mean(final_overlap_data(0.75*n_datapoints:n_datapoints));
        
        %{
        final_overlap_data_b = smooth(final_overlap_data_b, 10);
        slope_data_b = gradient(final_overlap_data_b, time_per_datapoint);
        final_velocity_data_b = smooth(slope_data_b, 100);
        max_velocities_b(i_seed, i_length) = max(abs(final_velocity_data_b)) * 1000;
        %}
    end
end

% Average over all seeds
for i_length=1:n_lengths
    avg_final_lengths(i_length) = mean(final_lengths(:, i_length));
    avg_max_velocities(i_length) = mean(max_velocities(:, i_length));
    % get standard devation of velocities
    err_max_velocities(i_length) = std(max_velocities(:, i_length));  
    %{
    avg_max_velocities_b(i_length) = mean(max_velocities_b(:, i_length));
    err_max_velocities_b(i_length) = std(max_velocities_b(:, i_length));
    %}
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

xticks([0,1,2,3]); %,6,8,10]);
%yticks([0,10,20,30,40,50]);
xlim([0 3]);
%ylim([0 55]);
axes = gca;
axes.XAxis.FontSize = 14;
axes.YAxis.FontSize = 14;
set(gcf, 'color', 'w');