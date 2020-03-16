clear all;
% READ FIXME BELOW
% Parameters from sim
mt_lengths = [1000, 500];
simName = 'test';
n_steps = 60000000;
n_datapoints = 10000;
delta_t = 0.0000025; 
unpin_time = 100;

n_sites_max = max(mt_lengths);
n_sites_min = min(mt_lengths);
n_mts = length(mt_lengths);
steps_per_datapoint = n_steps / n_datapoints;

xlinkSpeciesID = 1;
motorSpeciesID = 2;
start_time = 0;
end_time = n_steps * delta_t;

loyolagreen = 1/255*[0,104,87];
loyolagray = 1/255*[200,200,200];

% File info
fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
mtCoordFileName = '%s_mt_coord.file';
occupancyFileName = '%s_occupancy.file';

% Create filenames and open files
mtFile = sprintf(fileDirectory, sprintf(mtCoordFileName, simName));
occupancyFile = sprintf(fileDirectory, sprintf(occupancyFileName, simName));

% Read in and reshape data structure
mt_coord_data_file = fopen(mtFile);
mt_coord_raw_data = fread(mt_coord_data_file, n_mts * n_datapoints, '*int');
fclose(mt_coord_data_file);
mt_coord_data = reshape(mt_coord_raw_data, n_mts, n_datapoints);

occupancy_data_file = fopen(occupancyFile);
occupancy_raw_data = fread(occupancy_data_file, n_mts * n_sites_max * n_datapoints, '*int');
fclose(occupancy_data_file);
occupancy_data = reshape(occupancy_raw_data, n_sites_max, n_mts, n_datapoints);

overlap_data = zeros([n_datapoints 1]);
motor_frac_overlap = zeros([n_datapoints 1]);
motor_frac_nonoverlap = zeros([n_datapoints 1]);
xlink_frac_overlap = zeros([n_datapoints 1]);
xlink_frac_nonoverlap = zeros([n_datapoints 1]);

for i_data = 1:1:n_datapoints
    mt_coords = mt_coord_data(:, i_data);
    if mt_coords(1) > mt_coords(2)
        overlap_start = mt_coords(1);
    else
        overlap_start = mt_coords(2);
    end
    if((mt_coords(1) + mt_lengths(1)) > (mt_coords(2) + mt_lengths(2)))
        overlap_end = mt_coords(2) + mt_lengths(2);
    else
        overlap_end = mt_coords(1) + mt_lengths(1);    
    end 
    
    overlap_length = overlap_end - overlap_start; 
    overlap_data(i_data) = overlap_length;
    % below: number of SITES occupied by motors/xlinks (2 per double-bound)
    n_motors_overlap = 0;
    n_motors_nonoverlap = zeros(n_mts, 1);
    n_xlinks_overlap = 0;
    n_xlinks_nonoverlap = zeros(n_mts, 1);
    for i_mt = 1:1:n_mts
        occupancy = occupancy_data(:, i_mt, i_data);  
        n_sites = mt_lengths(i_mt);
        mt_coord = mt_coords(i_mt);
        for i_site = 1:1:n_sites
            site_coord = i_site + mt_coord;
            if occupancy(i_site) == motorSpeciesID
                if site_coord >= overlap_start && site_coord <= overlap_end
                    n_motors_overlap = n_motors_overlap + 1;
                else
                    n_motors_nonoverlap(i_mt) = n_motors_nonoverlap(i_mt) + 1;
                end
            elseif occupancy(i_site) == xlinkSpeciesID
                if site_coord >= overlap_start && site_coord <= overlap_end
                    n_xlinks_overlap = n_xlinks_overlap + 1;
                else
                    n_xlinks_nonoverlap(i_mt) = n_xlinks_nonoverlap(i_mt) + 1;
                end
            end
        end
    end
    motor_frac = n_motors_overlap / (n_mts * overlap_length);
    xlink_frac = n_xlinks_overlap / (n_mts * overlap_length);
    motor_frac_outside = 0.0;
    xlink_frac_outside = 0.0;
    for i_mt=1:n_mts
        delta = abs(mt_lengths(i_mt) - overlap_length);
        if delta ~= 0
            motor_frac_outside = motor_frac_outside + (n_motors_nonoverlap(i_mt) / delta);
            xlink_frac_outside = xlink_frac_outside + (n_xlinks_nonoverlap(i_mt) / delta);        
        end
    end
    motor_frac_overlap(i_data) = motor_frac;
    motor_frac_nonoverlap(i_data) = motor_frac_outside;
    xlink_frac_overlap(i_data) = xlink_frac;
    xlink_frac_nonoverlap(i_data) = xlink_frac_outside;
end

overlap_data = smooth(overlap_data) * 0.008;
motor_frac_overlap = smooth(motor_frac_overlap);
motor_frac_nonoverlap = smooth(motor_frac_nonoverlap);
xlink_frac_overlap = smooth(xlink_frac_overlap);
xlink_frac_nonoverlap = smooth(xlink_frac_nonoverlap);

%xlink_overlap_avg = mean(xlink_frac_overlap(0.75*n_datapoints, :))
%xlink_nonoverlap_avg = mean(xlink_frac_nonoverlap(0.75*n_datapoints, :))

fig1 = figure();
set(fig1, 'Position', [50, 50, 2.5*480, 2.5*300])

% Plot overlap length data
subplot(3, 1, 1)
plot(linspace(start_time, end_time, n_datapoints), overlap_data, ...
    'LineWidth', 2);
line([unpin_time unpin_time], ylim, 'Linestyle', '--', 'Color', 'red');
title('Overlap length over time');
ylabel('Overlap length (um)');
%legend('Overall overlap', 'location', 'northeastoutside');
xlim([0 end_time]);
%ylim([0 max(overlap_data)]);

% Plot motor occupancy data
subplot(3, 1, 2)
plot(linspace(start_time, end_time, n_datapoints), motor_frac_overlap, ...
    'LineWidth', 2);
hold on
plot(linspace(start_time, end_time, n_datapoints), motor_frac_nonoverlap, ...
    'LineWidth', 2);
title('Motor occupancy');
ylabel({'Fraction of', 'sites occupied'});
legend('In overlap', 'Outside overlap', 'location', 'northwest');
%ylim([0 0.15]);
xlim([0 end_time]);

% Plot xlink occupancy data
subplot(3, 1, 3)
plot(linspace(start_time, end_time, n_datapoints), xlink_frac_overlap, ...
    'LineWidth', 2);
hold on
plot(linspace(start_time, end_time, n_datapoints), xlink_frac_nonoverlap, ...
    'LineWidth', 2);
line([unpin_time unpin_time], ylim, 'Linestyle', '--', 'Color', 'red');
%line(xlim, [xlink_overlap_avg xlink_overlap_avg], 'Linestyle', '--', 'Color', 'black');
%line(xlim, [xlink_nonoverlap_avg xlink_nonoverlap_avg], 'Linestyle', '--', 'Color','black');
title('Crosslinker occupancy')
xlabel('Time (s)');
ylabel('Fraction of sites occupied');
legend('In overlap', 'Outside overlap', 'location', 'northwest');
xlim([0 end_time]);