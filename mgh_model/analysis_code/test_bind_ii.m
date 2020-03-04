clear all;
baseName = 'test_overlap4b';
mt_lengths = [101,101];
mt_coords = [0, 0];
mt_endpoints = mt_coords + mt_lengths;
n_datapoints = 100000;
speciesID = 1;
max_length = max(mt_lengths);
n_mts = length(mt_lengths);

k_on = 0.000238; % 1/(nM*s)
c_bind = 4500;  % nM
k_off_ii = 14.3; % 1/s
k_d = k_off_ii / k_on;
avg_occu = c_bind / (c_bind + k_d);

kbT = 4.114; % pN * nm
f_occu = @(x) c_bind / (c_bind + k_d*exp(delta_u(x)/kbT));


fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
fileStruct = '%s_occupancy.file';
data_file = fopen(sprintf(fileDirectory, sprintf(fileStruct, baseName)));
raw_data = fread(data_file, [max_length n_mts*n_datapoints], '*int');
fclose(data_file);
raw_data(raw_data ~= speciesID) = 0;
raw_data(raw_data == speciesID) = 1;

mt_one = zeros(1, mt_lengths(1));
mt_two = zeros(1, mt_lengths(2));

% Avg occupancy data for each MT over all datapoints
for i_data=1:2:((2*n_datapoints)-1)
    mt_one = mt_one + double(raw_data(1:mt_lengths(1), i_data)')./n_datapoints;
    mt_two = mt_two + double(raw_data(1:mt_lengths(2), i_data + 1)')./n_datapoints;
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
for i_mt=1:n_mts
    leftover_length = mt_lengths(i_mt) - overlap_length;
    for i_site=1:sys_size
        % If inside overlap, add it to overlap avg
        if i_site >= overlap_start && i_site <= overlap_end
            avg_occu_inside = avg_occu_inside + final_data(i_mt, i_site)/(2*overlap_length);
        elseif leftover_length ~= 0
            avg_occu_outside = avg_occu_outside + final_data(i_mt, i_site)/(2*leftover_length);
        end
    end
end

%%plot fig%%
fig1 = figure();
set(fig1,'Position', [50, 50, 2.5*480, 2.5*300])
plot(linspace(0, sys_size, sys_size), final_data(2, :), 'LineWidth', 2);
hold on
%plot([0 sys_size], [avg_occu avg_occu]);
fplot(@(x) f_occu(x - 50.5));


%%style stuff%%
ylim([0 avg_occu]);
xlim([0 sys_size]);
xlabel('Site coordinate (8 nm each)');

function u = delta_u(x)
    k_spring = 0.453; % pN / nm
    site_size = 8; % nm
    dr = sqrt(32^2 + (x*site_size)^2) - 32;
    u = 0.5*k_spring*dr*dr;
end
