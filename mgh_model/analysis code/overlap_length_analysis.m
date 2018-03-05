clear all
n_mts = 2;
n_sites = 625;
length = 625 * 8;
n_datapoints = 100000;
starting_point = 0000;
delta_t = 0.00001;
end_time = n_datapoints * delta_t * 200;
start_time = starting_point * delta_t * 200;

final_data = zeros([n_datapoints 1]);

fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
fileName = 'YES_MTcoord.file';

data_file = fopen(sprintf(fileDirectory, fileName));
raw_data = fread(data_file, [n_mts, n_datapoints], '*int');
fclose(data_file);

for i=starting_point+1:1:n_datapoints;
    mt_coord_one = raw_data(1, i);
    mt_coord_two = raw_data(2, i);
    delta = (mt_coord_two - mt_coord_one) * 8;
    overlap_length = length - delta;
    final_data(i, 1) = overlap_length; 
end

fig1 = figure();
set(fig1, 'Position', [50, 50, 2.5*480, 2.5*300])
plot(linspace(0, end_time, n_datapoints), final_data);

grid on
grid minor
title(...
        sprintf('Overlap Length Over Time (%g microns or %d sites in length)', ...
        n_sites * 8 / 1000, n_sites));
xlabel('Time (s)');
ylabel('Overlap length (nm)');
xlim([start_time end_time])