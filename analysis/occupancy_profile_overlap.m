clear all;
baseName = 'test_overlap4b';
mt_lengths = [101, 101];
mt_coords = [0, 0];
mt_endpoints = mt_coords + mt_lengths;
n_datapoints = 100000;
speciesID = 1;
max_length = max(mt_lengths);
n_mts = length(mt_lengths);

fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
fileStruct = '%s_occupancy.file';
data_file = fopen(sprintf(fileDirectory, sprintf(fileStruct, baseName)));
raw_data = fread(data_file, [max_length n_mts * n_datapoints], '*int');
fclose(data_file);
raw_data(raw_data ~= speciesID) = 0;
raw_data(raw_data == speciesID) = 1;

mt_one = zeros(1, mt_lengths(1));
mt_two = zeros(1, mt_lengths(2));

% Avg occupancy data for each MT over all datapoints
for i_data = 1:2:((2 * n_datapoints) - 1)
    mt_one = mt_one + double(raw_data(1:mt_lengths(1), i_data)') ./ n_datapoints;
    mt_two = mt_two + double(raw_data(1:mt_lengths(2), i_data + 1)') ./ n_datapoints;
end

% To account for MTs that are offset by some amount, create a larger
% 'system' that can hold all the data
sys_start = min(mt_coords);
sys_end = max(mt_endpoints);
sys_size = sys_end - sys_start;
final_data = zeros(n_mts, sys_size);

adjusted_coords = mt_coords - sys_start + 1;
adjusted_endpoints = mt_endpoints - sys_start;

final_data(1, adjusted_coords(1):adjusted_endpoints(1)) = mt_one;
final_data(2, adjusted_coords(2):adjusted_endpoints(2)) = mt_two;

% Find overlap region to calculate average occupances inside & outside it
overlap_start = max(adjusted_coords);
overlap_end = min(adjusted_endpoints);
overlap_length = overlap_end - overlap_start;

avg_occu_inside = 0.0;
avg_occu_outside = 0.0;

for i_mt = 1:n_mts
    leftover_length = mt_lengths(i_mt) - overlap_length;

    for i_site = 1:sys_size
        % If inside overlap, add it to overlap avg
        if i_site >= overlap_start && i_site <= overlap_end
            avg_occu_inside = avg_occu_inside + final_data(i_mt, i_site) / (2 * overlap_length);
        elseif leftover_length ~= 0
            avg_occu_outside = avg_occu_outside + final_data(i_mt, i_site) / (2 * leftover_length);
        end

    end

end

inside_occu_report = sprintf('\nAverage occupancy inside overlap: %g', avg_occu_inside);
%disp(inside_occu_report);
outside_occu_report = sprintf('Average occupancy outside overlap: %g', avg_occu_outside);
%disp(outside_occu_report);
bias = avg_occu_inside / avg_occu_outside;
bias_report = sprintf('Overlap binding bias: %g\n', bias);
%disp(bias_report);

%%plot fig%%
fig1 = figure();
set(fig1, 'Position', [50, 50, 2.5 * 480, 2.5 * 300])
plot(linspace(0, sys_size, sys_size), final_data(2, :), 'LineWidth', 2);
hold on
%plot([0 sys_size], [avg_occu_inside avg_occu_inside], '--r', 'Linewidth', 1.5);
%plot([0 sys_size], [avg_occu_outside avg_occu_outside], '--m', 'Linewidth', 1.5);

dim = [0.14 0.78 .1 .1];
str = {sprintf('Average overlap occupancy = %#.3g', avg_occu_inside), ...
        sprintf('Average non-overlap occupancy = %#.3g', avg_occu_outside), ...
        sprintf('Overlap binding bias = %#.3g', bias)};
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on');

%%style stuff%%
ylim([0 1]);
xlim([0 sys_size]);
title('c_{bulk} = 0.5 nM; k_{on} = 0.005 nM^{-1}s^{-1}; c_{eff} = 10 nM; L_{MT} = 2 microns; L_{overlap} = 1 micron');
xlabel('Site coordinate (8 nm each)');
ylabel('Fraction of the time occupied');
legend({'First MT occupancy', 'Second MT occupancy', 'Avg. overlap occupancy', 'Avg. non-overlap occupancy'}, ...
    'Location', 'northeastoutside');
