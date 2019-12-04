clear all;
% Effective binding concentration for tethered motors/xlinks
c_eff_teth = 4500;
% Effective binding concentration for second head of xlinks
c_eff_bind = 4500;
% Parameter to be scanned over for multiple scaling plots
scan_param = [1]; % NULL
legendLabel = cellstr(num2str(scan_param','k off ratio = %i'));
% Length of top, mobile MT in n_sites
top_lengths = [75,150,225,300,375,450,525,600];
seeds = [1]; %,2,3,4,5,6,7,8,9,11,10,12];
base_name = "2019_12_02_slideScan/slide_scan_%i_%i";
fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
site_size = 0.008;

final_length = zeros(length(scan_param), length(top_lengths), length(seeds));
for i_scan=1:length(scan_param)
    for i_length=1:length(top_lengths)
        for i_seed=1:length(seeds)
            sim_name = sprintf(base_name, top_lengths(i_length), seeds(i_seed))
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
            n_datapoints = str2double(values{contains(params, "n_datapoints")});
            % Open mt coordinate file
            mt_coords_file = fopen(sprintf(fileDirectory, sprintf('%s_mt_coord.file', sim_name)));
            mt_coord_data = fread(mt_coords_file, [n_mts, n_datapoints], 'double');
            fclose(mt_coords_file);
            
            overlap_length = zeros(n_datapoints, 1);
            for i_data=1:n_datapoints
                mt_coords = mt_coord_data(:, i_data)';
                mt_endpoints = mt_coords + mt_lengths;
                overlap_start = max(mt_coords);
                overlap_end = min(mt_endpoints);
                overlap_length(i_data) = (overlap_end - overlap_start) * site_size;
            end
            
            avg_length = mean(overlap_length(0.9*n_datapoints:n_datapoints))
            
            final_length(i_scan, i_length, i_seed) = avg_length;
        end
    end
end

final_data = final_length(:, :, 1);

fig1 = figure();
set(fig1, 'Position', [50, 50, 640, 420]);
plot(top_lengths*site_size, final_data, '-o', 'LineWidth', 2);
title(sprintf('Sliding with c eff teth = %i, c eff bind = %i', c_eff_teth, c_eff_bind), ...
    'FontSize', 14);
ylabel('Final overlap lenth (microns)', 'FontName', 'Garuda','FontSize', 12);
xlabel('Initial overlap length (microns)', 'FontName', 'Garuda','FontSize', 12)
axes = gca;
axes.XAxis.FontSize = 14;
axes.YAxis.FontSize = 14;
if length(scan_param) > 1
    legend(legendLabel, 'Location', 'best');
end