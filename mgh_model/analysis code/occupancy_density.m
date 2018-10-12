clear all

% Parameters from sim
n_steps = 2000000;
n_sites = 500;
n_mts = 2;
n_datapoints = 100000;
delta_t = 0.0005; 
xlinkSpeciesID = 1;
motorSpeciesID = 2;
start_time = 0;
end_time = n_steps * delta_t;

% File info
simName = 'slide_testlong';
fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
mtFileName = '%s_mt_coord.file';
occupancyFileName = '%s_occupancy.file';

% Create filenames and open files
mtFile = sprintf(fileDirectory, sprintf(mtFileName, simName));
occupancyFile = sprintf(fileDirectory, sprintf(occupancyFileName, simName));

% Read in and reshape data structure
mt_data_file = fopen(mtFile);
mt_raw_data = fread(mt_data_file, [n_mts * n_datapoints], '*double');
fclose(mt_data_file);
mt_data = reshape(mt_raw_data, n_mts, n_datapoints);

occupancy_data_file = fopen(occupancyFile);
occupancy_raw_data = fread(occupancy_data_file, [n_mts * n_sites * n_datapoints], '*int');
fclose(occupancy_data_file);
occupancy_data = reshape(occupancy_raw_data, n_sites, n_mts, n_datapoints);

overlap_data = zeros([n_datapoints 1]);
motor_frac_overlap = zeros([n_datapoints 1]);
motor_frac_nonoverlap = zeros([n_datapoints 1]);
xlink_frac_overlap = zeros([n_datapoints 1]);
xlink_frac_nonoverlap = zeros([n_datapoints 1]);

for i_data = 1:1:n_datapoints
    mt_coords = mt_data(:, i_data);
    if mt_coords(1) > mt_coords(2)
        overlap_start = mt_coords(1);
        overlap_end = mt_coords(2) + n_sites;
    else
        overlap_start = mt_coords(2);
        overlap_end = mt_coords(1) + n_sites;
    end
    overlap_length = overlap_end - overlap_start; 
    overlap_data(i_data) = overlap_length;
    
    % below: number of SITES occupied by motors/xlinks (2 per double-bound)
    n_motors_overlap = 0;
    n_motors_nonoverlap = 0;
    n_xlinks_overlap = 0;
    n_xlinks_nonoverlap = 0;
    
    for i_mt = 1:1:n_mts
        occupancy = occupancy_data(:, i_mt, i_data);
        mt_coord = mt_coords(i_mt);
        for i_site = 1:1:n_sites
            site_coord = i_site + mt_coord;
            if occupancy(i_site) == motorSpeciesID
                if site_coord >= overlap_start && site_coord <= overlap_end
                    n_motors_overlap = n_motors_overlap + 1;
                else
                    n_motors_nonoverlap = n_motors_nonoverlap + 1;
                end
            elseif occupancy(i_site) == xlinkSpeciesID
                if site_coord >= overlap_start && site_coord <= overlap_end
                    n_xlinks_overlap = n_xlinks_overlap + 1;
                else
                    n_xlinks_nonoverlap = n_xlinks_nonoverlap + 1;
                end
            end
        end
    end
    motor_frac = n_motors_overlap / (n_mts * overlap_length);
    xlink_frac = n_xlinks_overlap / (n_mts * overlap_length);
    if overlap_length ~= n_sites
        motor_frac_outside = n_motors_nonoverlap / (n_mts*(n_sites - overlap_length));
        xlink_frac_outside = n_xlinks_nonoverlap / (n_mts*(n_sites - overlap_length));
    else
        motor_frac_outside = 0; 
        xlink_frac_outside = 0;
    end

    motor_frac_overlap(i_data) = motor_frac;
    motor_frac_nonoverlap(i_data) = motor_frac_outside;
    xlink_frac_overlap(i_data) = xlink_frac;
    xlink_frac_nonoverlap(i_data) = xlink_frac_outside;
end

fig1 = figure();
set(fig1, 'Position', [50, 50, 2.5*480, 2.5*300])

% Plot overlap length data
subplot(3, 1, 1)
plot(linspace(start_time, end_time, n_datapoints), overlap_data * 0.008, ...
    'LineWidth', 2);
title('Overlap length in microns');
legend('Overall overlap', 'location', 'northeastoutside');

% Plot motor occupancy data
subplot(3, 1, 2)
plot(linspace(start_time, end_time, n_datapoints), motor_frac_overlap, ...
    'LineWidth', 2);
hold on
plot(linspace(start_time, end_time, n_datapoints), motor_frac_nonoverlap, ...
    'LineWidth', 2);
title('Motor occupancy');
legend('In overlap', 'Outside overlap', 'location', 'northeastoutside');

% Plot xlink occupancy data
subplot(3, 1, 3)
plot(linspace(start_time, end_time, n_datapoints), xlink_frac_overlap, ...
    'LineWidth', 2);
hold on
plot(linspace(start_time, end_time, n_datapoints), xlink_frac_nonoverlap, ...
    'LineWidth', 2);
title('Crosslinker occupancy')
legend('In overlap', 'Outside overlap', 'location', 'northeastoutside');
