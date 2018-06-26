clear all
n_datapoints = 100000;
species_ID = 2;

k_off_ratio = 2500;
failstep_rates = [0.005, 0.05, 0.5, 5, 25, 50];
k_off_pseudos = [0.5, 1, 5, 10, 50, 100];
mt_lengths = [250, 500, 750, 1000, 1250];

failstep_axis_ticks = cellstr(num2str(failstep_rates', '%g'));
k_off_legend = cellstr(num2str(k_off_pseudos', 'k-off pseudo = %-g'));

% Set color order so lines have consistent colors for up to 8 k_off_pseudo values
ColorArray = {'orange', 'cyan', 'blue', 'green', 'yellow', ... 
    'dark grey','magenta', 'red'};
RGBArray = [0.91 0.41 0.17; 0.301 0.745 0.933; 0 0 1; 0 0.5 0; ...
    0.929 0.6940 0.125; 0.3 0.3 0.3; 1 0 1; 1 0 0];

fileDirectory = '/media/shane/Shane''s External HDD (1 TB)/Parameter Scan 1/%s';
%fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
nameStruct = '%i_%#.3f_%#.2f_%i_occupancy.file';


endtag_data = zeros([length(k_off_pseudos) length(mt_lengths) length(failstep_rates) 1]);

% Run through different failstep rates
for i_fsr=1:1:1 %length(failstep_rates)
    failstep_rate = failstep_rates(i_fsr);
    % For each failstep rate, run through different k_off_pseudo values
    for j_kop=1:1:length(k_off_pseudos)
        k_off_pseudo = k_off_pseudos(j_kop);
        % Scan through respective MT lengths for each param set
        for k_mtl = 1:1:length(mt_lengths)
            mt_length = mt_lengths(k_mtl);
            
            temp_mt = zeros([mt_length 1]);
            final_mt = zeros([mt_length 1]);
            
            fileName = sprintf(nameStruct, k_off_ratio, failstep_rate, ...
                k_off_pseudo, mt_length)
            
            data_file = fopen(sprintf(fileDirectory, fileName));
            raw_data = fread(data_file, [mt_length, n_datapoints], '*int');
            fclose(data_file);
            
            raw_data((raw_data ~= species_ID) | (raw_data == 0)) = 0;
            raw_data((raw_data == species_ID) | (raw_data ~= 0)) = 1;
            shape = InsertShape(r);
            % Read in and average occupancy data over all datapoints
            for i=1:1:n_datapoints
                temp_mt(:) = raw_data(:,i);
                final_mt(:) = final_mt(:) + double(temp_mt(:)./n_datapoints);
            end
            
            % Find endtag length based on averaged (final) MT occupancy data
            endtag_site = mt_length;
            for i=1:1:mt_length
                site_occupancy = final_mt(i, 1);
                if(site_occupancy > 0.5)
                    endtag_site =  i;
                    break;
                end
            end
            endtag_length = (mt_length - endtag_site) * 0.008
            endtag_data(j_kop, k_mtl, i_fsr) = endtag_length
        end
    end
end

% Plot data

% x-axis always corresponds to MT length values
x_data = mt_lengths / 125;
% y-axis corresponds to different failstep_rate values
for i_y=1:1:length(failstep_rates)
    % Make sure each failstep_rate "slice" has a constant y-value
    y = failstep_rates(i_y);
    y_data(1:length(mt_lengths)) = i_y;
    % z-axis will be avg occupancy 
    for i_z=1:1:length(k_off_pseudos)
        z_data = endtag_data(i_z, :, i_y);
        plot3(x_data, y_data, z_data, 'Color', RGBArray(i_z, :), ...
            'LineWidth', 1.5);
        hold on
    end
end

% Style stuff
view(-80, 15)
set(gcf, 'Position', [50 50 1080 800]);
grid on 
grid minor
legend(k_off_legend, 'Location', 'northeastoutside');
title(sprintf('k-off ratio = %g (stalled motors unbind %gx less)', ... 
    k_off_ratio, k_off_ratio));
xlabel("Microtubule length (microns)");
%xlabel_units = get(get(gca, 'XLabel'), 'Units');
%set(get(gca, 'XLabel'), 'Units', 'normalized');
xlabel_pos = get(get(gca, 'XLabel'), 'Position');
set(get(gca, 'XLabel'), 'Position', xlabel_pos + [0 -1.3 0.15]);
%set(get(gca, 'XLabel'), 'Units', xlabel_units);
ylabel("Failstep rate (per second for stalled motors)"); 
yticks(1:1:length(failstep_rates));
yticklabels(failstep_axis_ticks);
zlabel("Length of endtag (microns)");