clear variables;

sim_name = 'xlFreeze'; % Raw sim name; do not include directory
output_movie_name = 'mov';

start_frame = 1;
frames_per_plot = 100; % in n_datapoints; number of timesteps per output plot
end_frame = -1;  % set to -1 to run until end of data
movie_duration = 30; % in seconds

endtag_region_size = 50;  % distance (in n_sites) from plus-end that defines end-tag region

% Load parameter structure
file_dir = '..';  % Default; only change if you move CyLaKS output files
params = load_parameters(sprintf('%s/%s', file_dir, sim_name));

% Set plotting helper variables
sid_site = 0;
sid_xlink = 1;
sid_motor = 2;
purple = [128 0 128] / 255; 
blue = [30 144 255] / 255;
if end_frame == -1
    end_frame = params.n_datapoints;
end
active_frames = end_frame - start_frame + 1;

% Initialize videowriter object
v = VideoWriter(output_movie_name);
v.FrameRate = (active_frames / frames_per_plot) / movie_duration;
open(v);
frame_box = [0 0 1445 300];

% Open figure and set to desired size (each frame must be this same size)
fig1 = figure;
set(fig1, 'Position', [15 15 1080 640]);

% Open data files
filament_filename = sprintf('%s/%s_filament_pos.file', file_dir, sim_name);
filament_pos = zeros(params.n_dims, 2, params.n_mts, params.n_datapoints);
filament_pos = load_data(filament_pos, filament_filename, '*double');

protein_filename = sprintf('%s/%s_protein_id.file', file_dir, sim_name);
protein_ids = zeros(params.max_sites, params.n_mts, params.n_datapoints) - 1;
protein_ids = load_data(protein_ids, protein_filename, '*int');

occupancy_filename = sprintf('%s/%s_occupancy.file', file_dir, sim_name);
occupancy = zeros(params.max_sites, params.n_mts, params.n_datapoints) - 1;
occupancy = load_data(occupancy, occupancy_filename, '*int');

partner_filename = sprintf('%s/%s_partner_index.file', file_dir, sim_name);
partner_indices = zeros(params.max_sites, params.n_mts, params.n_datapoints) - 1;
partner_indices = load_data(partner_indices, partner_filename, '*int');

xlink_filename = sprintf('%s/%s_xlink_force.file', file_dir, sim_name);
xlink_data = zeros(params.n_dims, params.n_datapoints);
xlink_data = load_data(xlink_data, xlink_filename, 'double');

n_sites_top = params.mt_lengths(2);
n_sites_bot = params.mt_lengths(1);

% Main data structures -- construct before loops to enhance speed
% (full data)
num_escape_left = zeros(active_frames, 1);
num_double_left = zeros(active_frames, 1);
num_double_bulk = zeros(active_frames, 1);
num_double_right = zeros(active_frames, 1);
num_escape_right = zeros(active_frames, 1);
num_single_bulk = zeros(active_frames, 1);
% (time-averaged data)
n_escape_left = zeros(active_frames/frames_per_plot, 1) - 100;
n_double_left = zeros(active_frames/frames_per_plot, 1) - 100;
n_double_bulk = zeros(active_frames/frames_per_plot, 1) - 100;
n_double_right = zeros(active_frames/frames_per_plot, 1) - 100;
n_escape_right = zeros(active_frames/frames_per_plot, 1) - 100;
n_single_bulk = zeros(active_frames/frames_per_plot, 1) - 100;
% Main loop -- each iteration will update plots and be a movie frame
for i_frame = start_frame : frames_per_plot : end_frame - frames_per_plot
    % Sub loop -- allows us to bin all data between movie frames
    n_single = 0;
    n_double = 0;
    %single_coords = zeros(2*params.max_sites, params.n_dims);
    %double_coords = zeros(params.max_sites, 2, params.n_dims);
    angles = zeros(params.max_sites, 1);
    lengths = zeros(params.max_sites, 1);
    neighb_distance = zeros(2*params.max_sites, 1);
    

    for i_subframe = 0 : 1 : frames_per_plot
        i_timestep = i_frame + i_subframe;
        % Extract microtubule tip coords from filament_pos data
        plus_pos_top = filament_pos(:, 1, 2, i_timestep);
        minus_pos_top = filament_pos(:, 2, 2, i_timestep);
        plus_pos_bot = filament_pos(:, 1, 1, i_timestep);
        minus_pos_bot = filament_pos(:, 2, 1, i_timestep);
        % Use tip coords to calculate overlap boundaries
        overlap_start = max(plus_pos_bot(1), minus_pos_top(1));
        left_endtag_end = overlap_start + endtag_region_size * 8.2;
        overlap_end = min(plus_pos_top(1), minus_pos_bot(1));
        right_endtag_start = overlap_end  - endtag_region_size * 8.2;
        % Get variables that allow us to go from site index to coordinates
        top_mt_vec = [plus_pos_top(1) - minus_pos_top(1), plus_pos_top(2) - minus_pos_top(2)];
        bot_mt_vec = [minus_pos_bot(1) - plus_pos_bot(1), minus_pos_bot(2) - plus_pos_bot(2)];
        % Scan over top (i_mt = 2) microtubule for singly-bound ONLY
        for i_site = 1 : 1 : n_sites_top
            if occupancy(i_site, 2, i_timestep) == sid_xlink
                % If partner index doesn't exist, it's singly bound
                if partner_indices(i_site, 2, i_timestep) == -1
                    pos_x = minus_pos_top(1) + (double(i_site - 1)/(n_sites_top - 1))*top_mt_vec(1);
                    pos_y = minus_pos_top(2) + (double(i_site - 1)/(n_sites_top - 1))*top_mt_vec(2);
                    %n_single = n_single + 1;
                    %single_coords(n_single, :) = [pos_x, pos_y];
                    if pos_x < overlap_start
                        num_escape_left(i_timestep) = num_escape_left(i_timestep) + 1;
                    elseif pos_x > overlap_end
                        num_escape_right(i_timestep) = num_escape_right(i_timestep) + 1;
                    else
                        num_single_bulk(i_timestep) = num_single_bulk(i_timestep) + 1;
                    end
                end
            end
        end
        % Scan over bottom (i_mt = 1) microtutuble
        for i_site = 1 : 1 : n_sites_bot
            if occupancy(i_site, 1, i_timestep) == sid_xlink
                pos_x = plus_pos_bot(1) + (double(i_site - 1)/(n_sites_bot - 1))*bot_mt_vec(1);
                pos_y = plus_pos_bot(2) + (double(i_site - 1)/(n_sites_bot - 1))*bot_mt_vec(2);
                % If partner index doesn't exist, it's singly bound
                if partner_indices(i_site, 1, i_timestep) == -1
                    %n_single = n_single + 1;
                    %single_coords(n_single, :) = [pos_x, pos_y];
                    if pos_x < overlap_start
                        num_escape_left(i_timestep) = num_escape_left(i_timestep) + 1;
                    elseif pos_x > overlap_end
                        num_escape_right(i_timestep) = num_escape_right(i_timestep) + 1;
                    else
                        num_single_bulk(i_timestep) = num_single_bulk(i_timestep) + 1;
                    end
                else
                    ii_site = partner_indices(i_site, 1, i_timestep);
                    endpos_x = minus_pos_top(1) + (double(ii_site)/(n_sites_top - 1))*top_mt_vec(1);
                    endpos_y = minus_pos_top(2) + (double(ii_site)/(n_sites_top - 1))*top_mt_vec(2);
                    len_x = endpos_x - pos_x;
                    len_y = endpos_y - pos_y;
                    len = sqrt(len_x^2 + len_y^2);
                    theta = acos(len_x/len) * 180 / pi;  % get angle and convert to degrees
                    % record statistics
                    n_double = n_double + 1;
                    %double_coords(n_double, 1, :) = [pos_x, pos_y];
                    %double_coords(n_double, 2, :) = [endpos_x, endpos_y];
                    angles(n_double) = theta;
                    lengths(n_double) = len;
                    if pos_x <= left_endtag_end || endpos_x <= left_endtag_end
                       num_double_left(i_timestep) = num_double_left(i_timestep) + 1;    
                    end
                    if pos_x >= right_endtag_start || endpos_x >= right_endtag_start
                       num_double_right(i_timestep) = num_double_right(i_timestep) + 1; 
                    end
                    if (pos_x > left_endtag_end && pos_x < right_endtag_start) && ...
                       (endpos_x > left_endtag_end && endpos_x < right_endtag_start) && ...
                       left_endtag_end < right_endtag_start
                       num_double_bulk(i_timestep) = num_double_bulk(i_timestep) + 1;
                    end
                    % find distance to nearest neighbor for each head
                    for i_mt = 1 : 1 : params.n_mts
                        neighb_found = false;
                        for delta = 1 : 1 : params.max_sites
                            if neighb_found
                                break;
                            end
                            for dir = -1 : 2 : 1
                                i_scan = i_site + (dir * delta);
                                if i_scan >= 1 && i_scan <= params.mt_lengths(i_mt)
                                    sid_neighb = occupancy(i_scan, i_mt, i_timestep);
                                    if sid_neighb == sid_xlink
                                        neighb_distance(2*n_double + i_mt - 1) = delta - 1;
                                        neighb_found = true;
                                        break;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    % Trim data
    angles = angles(1:n_double);
    lengths = lengths(1:n_double);
    neighb_distance = neighb_distance(1:2*n_double);
    
    % Reset MT positions to index of frame we are plotting
    plus_pos_top = filament_pos(:, 1, 2, i_frame);
    minus_pos_top = filament_pos(:, 2, 2, i_frame);
    plus_pos_bot = filament_pos(:, 1, 1, i_frame);
    minus_pos_bot = filament_pos(:, 2, 1, i_frame);
    overlap_start = max(plus_pos_bot(1), minus_pos_top(1));
    overlap_end = min(plus_pos_top(1), minus_pos_bot(1));
    
    % Clear figure so that it only displays figures from current frame
    clf;
    % Plot occupancy profile
    subplot(2, 4, 1)
    occupancy(occupancy ~= sid_xlink) = 0;
    occupancy(occupancy == sid_xlink) = 1;
    avg_top = zeros([params.mt_lengths(2) 1]);
    avg_bot = zeros([params.mt_lengths(1) 1]);
    for i_window = 0 : 1 : frames_per_plot
        for i_site = 1 : params.mt_lengths(1)
            avg_top(i_site) = avg_top(i_site) + double(occupancy(i_site, 2, i_frame + i_window))./ frames_per_plot;
        end
        for i_site = 1 : params.mt_lengths(2)
            avg_bot(i_site) = avg_bot(i_site) + double(occupancy(i_site, 1, i_frame + i_window)) ./ frames_per_plot;
        end
    end
    x_begin = minus_pos_top(1);
    x_end = minus_pos_bot(1);
    span = x_end - x_begin;
    span_sites = int32(span / (params.site_size * 1000));
    i_offset = span_sites - params.mt_lengths(2);
    if i_offset < 1
        i_offset = 1;
    end
    avg_tot = zeros([span_sites 1]);
    avg_tot(1:params.mt_lengths(1)) = avg_top;
    avg_tot(i_offset:i_offset+params.mt_lengths(2) - 1) = avg_tot(i_offset:i_offset+params.mt_lengths(2) - 1 ) + avg_bot;
    plot(avg_tot, 'LineWidth', 1.5);
    xlim([0 span_sites]);
    ylim([-0.1 2]);
    ylabel("Fractional occupancy");
    xlabel("Site number");
    
    
    
    % Plot simplified movie snapshot
    subplot(2, 4, 2)
    hold all
    if overlap_start < overlap_end
        xlim([(overlap_start - 200) (overlap_end + 200)]);
    else
        xlim([(overlap_end - 200) (overlap_start + 200)]);
    end
    ylim([(plus_pos_bot(2) - 2.5) (plus_pos_top(2) + 2.5)]);
    line([plus_pos_top(1), minus_pos_top(1)],[plus_pos_top(2), minus_pos_top(2)], ...
        'LineWidth', 4, 'Color', [0.7 0.7 0.7]);
    line([plus_pos_bot(1), minus_pos_bot(1)],[plus_pos_bot(2), minus_pos_bot(2)], ...
        'LineWidth', 4, 'Color', [0.7 0.7 0.7]);
    for i_site = 1 : 1 : params.max_sites
        id = protein_ids(i_site, 1, i_frame);
        sid = occupancy(i_site, 1, i_frame);
        if sid == sid_xlink
            n_sites = params.mt_lengths(1);
            bot_mt_vec = [minus_pos_bot(1) - plus_pos_bot(1), minus_pos_bot(2) - plus_pos_bot(2)];
            pos_x = plus_pos_bot(1) + (double(i_site-1)/(n_sites-1))*bot_mt_vec(1);
            pos_y = plus_pos_bot(2) + (double(i_site-1)/(n_sites-1))*bot_mt_vec(2);
            if partner_indices(i_site, 1, i_frame) ~= -1
                ii_site = partner_indices(i_site, 1, i_frame);
                nn_sites = params.mt_lengths(2);
                top_mt_vec = [plus_pos_top(1) - minus_pos_top(1), plus_pos_top(2) - minus_pos_top(2)];
                endpos_x = minus_pos_top(1) + (double(ii_site)/(nn_sites-1))*top_mt_vec(1);
                endpos_y = minus_pos_top(2) + (double(ii_site)/(nn_sites-1))*top_mt_vec(2);
                line([pos_x, endpos_x],[pos_y, endpos_y], ...
                    'LineWidth', 1, 'Color', purple);
            end
        end
    end
    xlabel("Position in x (nm)");
    ylabel("Position in y (nm)");
    
    
    
    
    % Plot 1-to-1 axis plot so we can see angles accurately
    subplot(2, 4, 3)
    hold all
    %xlim([(overlap_start - 10) (overlap_start + 40)]);
    xlim([(overlap_end - 40) (overlap_end + 10)]);
    ylim([(plus_pos_bot(2) - 10) (plus_pos_bot(2) + 40)]);
    line([plus_pos_top(1), minus_pos_top(1)],[plus_pos_top(2), minus_pos_top(2)], ...
        'LineWidth', 4, 'Color', [0.7 0.7 0.7]);
    line([plus_pos_bot(1), minus_pos_bot(1)],[plus_pos_bot(2), minus_pos_bot(2)], ...
        'LineWidth', 4, 'Color', [0.7 0.7 0.7]);
    for i_site = 1 : 1 : params.max_sites
        id = protein_ids(i_site, 1, i_frame);
        sid = occupancy(i_site, 1, i_frame);
        if sid == sid_xlink
            n_sites = params.mt_lengths(1);
            bot_mt_vec = [minus_pos_bot(1) - plus_pos_bot(1), minus_pos_bot(2) - plus_pos_bot(2)];
            pos_x = plus_pos_bot(1) + (double(i_site-1)/(n_sites-1))*bot_mt_vec(1);
            pos_y = plus_pos_bot(2) + (double(i_site-1)/(n_sites-1))*bot_mt_vec(2);
            if partner_indices(i_site, 1, i_frame) ~= -1
                ii_site = partner_indices(i_site, 1, i_frame);
                nn_sites = params.mt_lengths(2);
                top_mt_vec = [plus_pos_top(1) - minus_pos_top(1), plus_pos_top(2) - minus_pos_top(2)];
                endpos_x = minus_pos_top(1) + (double(ii_site)/(nn_sites-1))*top_mt_vec(1);
                endpos_y = minus_pos_top(2) + (double(ii_site)/(nn_sites-1))*top_mt_vec(2);
                line([pos_x, endpos_x],[pos_y, endpos_y], ...
                    'LineWidth', 2.5, 'Color', purple);
            end
        end
    end
    xlabel("Position in x (nm)");
    ylabel("Position in y (nm)");
    
    
    
    
    % delta N plot here
    subplot(2, 4, 4)
    x_axis = [start_frame:frames_per_plot:end_frame-1];
    hold all
    
    n_escape_left(floor(i_frame/frames_per_plot) + 1) = mean(num_escape_left(i_frame:i_frame+frames_per_plot));
    n_double_left(floor(i_frame/frames_per_plot) + 1) = mean(num_double_left(i_frame:i_frame+frames_per_plot));
    n_double_bulk(floor(i_frame/frames_per_plot) + 1) = mean(num_double_bulk(i_frame:i_frame+frames_per_plot));
    n_double_right(floor(i_frame/frames_per_plot) + 1) = mean(num_double_right(i_frame:i_frame+frames_per_plot));
    n_escape_right(floor(i_frame/frames_per_plot) + 1) = mean(num_escape_right(i_frame:i_frame+frames_per_plot));
    n_single_bulk(floor(i_frame/frames_per_plot) + 1) = mean(num_single_bulk(i_frame:i_frame+frames_per_plot));
    plot(x_axis, n_double_left, '<', 'LineWidth', 1.5, 'Color', [0, 150, 255] / 255);
    plot(x_axis, n_double_bulk, 'sq', 'LineWidth', 1.5, 'Color', [128,  0, 128] / 255);
    plot(x_axis, n_double_right, '>', 'LineWidth', 1.5, 'Color', [210, 4, 45] / 255);
    plot(x_axis, n_escape_left, 'o', 'LineWidth', 1.5, 'Color', [167, 199, 231] / 255);
    plot(x_axis, n_single_bulk, 'o', 'LineWidth', 1.5, 'Color', [195, 177, 225] / 255);
    plot(x_axis, n_escape_right, 'o', 'LineWidth', 1.5, 'Color', [250, 160, 160] / 255);
    legend(["2HB (L)", "2HB (C)", "2HB (R)", "1HB (L)", "1HB (C)","1HB (R)"], 'location', 'best');
    ylim([0 inf]);
    xlim([start_frame end_frame]); 
    ylabel("Count (global)");
    xlabel("Timestep number");
    
    
    
    % Plot histogram of xlink angles
    subplot(2, 4, 5)
    n_bins = ceil(sqrt(length(angles)) / 4.0);
    histfit(angles, n_bins, 'kernel');
    xlim([35 145]);
    xticks([40 90 140]);
    ylabel("Count");
    xlabel("Crosslinker angle (deg)");
    
  
    
    
    % Plot histogram of xlink lengths
    subplot(2, 4, 6)
    n_bins = round(sqrt(length(lengths)) / 4.0);
    histfit(lengths, n_bins, 'kernel');
    xlim([15 50]);
    xticks([20 32 45]);
    ylabel("Count");
    xlabel("Crosslinker length (nm)");
    
    
    
    % Plot histogram of xlink neighbor distance
    subplot(2, 4, 7)
    n_bins = max(neighb_distance);
    h = histfit(neighb_distance, n_bins, 'exponential');
    xlim([-1 15]);
    xticks([0 4 8]);
    ylabel("Count");
    xlabel("Distance to neighbor (n sites)");
    
    
    
    % n_bound histogram here
    subplot(2, 4, 8)
    bar([num_escape_left(i_frame) num_double_left(i_frame) num_double_bulk(i_frame) ...
        num_double_right(i_frame) num_escape_right(i_frame) num_single_bulk(i_frame)]');
    ylabel("Count");
    xlabel("Populations");
    xticklabels(["1HB (L)", "2HB (L)", "2HB (C)", "2HB (R)", "1HB (R)", "1HB (C)"]);
    xtickangle(45);
    
    
    time = (i_frame - start_frame) * params.time_per_datapoint;
    str = sprintf('Time: %#.2f seconds', time);
    dim = [0.15 0.69 .3 .3];
    annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on');
    str = sprintf('Force: %#.2f pN', xlink_data(1, i_frame));
    dim = [0.35 0.69 0.3 0.3];
    annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on');
    drawnow();
    writeVideo(v, getframe(gcf));
    
    if overlap_start > overlap_end
        break
    end
end
close(v);