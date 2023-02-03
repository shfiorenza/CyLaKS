clear variables;

sim_name = 'testino2'; % Raw sim name; do not include directory
output_movie_name = 'mov';

start_frame = 1;
end_frame = -1;  % set to -1 to run until end of data

frames_per_plot = 100; % in n_datapoints
movie_duration = 30; % in seconds

% Load parameter structure
file_dir = '..';  % Default; only change if you move CyLaKS output files
params = load_parameters(sprintf('%s/%s', file_dir, sim_name));

% Set plotting helper variables
sid_site = 0;
sid_xlink = 1;
sid_motor = 2;
purple = [128 0 128] / 255;
if end_frame == -1
    end_frame = params.n_datapoints; 
end
active_frames = end_frame - start_frame;

% Initialize videowriter object
v = VideoWriter(output_movie_name);
v.FrameRate = (active_frames / frames_per_plot) / movie_duration;
open(v);
frame_box = [0 0 1445 300];

% Open figure and set to desired size (each frame must be this same size)
fig1 = figure;
set(fig1, 'Position', [15 15 1500 750]);

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

stats_double = zeros(1, active_frames) - 100;
stats_single_in = zeros(1, active_frames) - 100;
stats_single_out = zeros(1, active_frames) - 100;

end_vid = false;
% Run through all datapoints; each one is a frame in our movie
for i_data = start_frame : frames_per_plot : end_frame - frames_per_plot
    
    if end_vid
       break; 
    end
          
    % Clear figure so that it only displays figures from current datapoint
    clf;
    
    % scan over microtubules and bin crosslinkers appropriately
    n_single_in = 0;
    n_single_out = 0;
    n_double = 0;
    % just doubly bound for now
    for i_window = 0 : 1 : frames_per_plot
    % determine overlap start and end coordinates 
    plus_pos_top = filament_pos(:, 1, 2, i_data + i_window);
    minus_pos_top = filament_pos(:, 2, 2, i_data + i_window);
    plus_pos_bot = filament_pos(:, 1, 1, i_data + i_window);
    minus_pos_bot = filament_pos(:, 2, 1, i_data + i_window);
    overlap_start = max(plus_pos_bot(1), minus_pos_top(1));
    overlap_end = min(plus_pos_top(1), minus_pos_bot(1));
    if overlap_start > overlap_end
        end_vid = true;
        break
    end
        
    for i_site = 1 : 1 : params.max_sites
        id = protein_ids(i_site, 1, i_data + i_window);
        sid = occupancy(i_site, 1, i_data + i_window);
        if sid == sid_xlink
            n_sites = params.mt_lengths(1);
            bot_mt_vec = [minus_pos_bot(1) - plus_pos_bot(1), minus_pos_bot(2) - plus_pos_bot(2)];
            pos_x = plus_pos_bot(1) + ((i_site-1)/(n_sites-1))*bot_mt_vec(1);
            pos_y = plus_pos_bot(2) + ((i_site-1)/(n_sites-1))*bot_mt_vec(2);
            if partner_indices(i_site, 1, i_data + i_window) ~= -1    
                ii_site = partner_indices(i_site, 1, i_data + i_window);
                nn_sites = params.mt_lengths(2);
                top_mt_vec = [plus_pos_top(1) - minus_pos_top(1), plus_pos_top(2) - minus_pos_top(2)];
                endpos_x = minus_pos_top(1) + (double(ii_site)/(nn_sites-1))*top_mt_vec(1);
                endpos_y = minus_pos_top(2) + (double(ii_site)/(nn_sites-1))*top_mt_vec(2);
                len_x = endpos_x - pos_x;
                len_y = endpos_y - pos_y;
                len = sqrt(len_x^2 + len_y^2);
                theta = acos(len_x/len) * 180 / pi;  % get angle and convert to degrees
                n_double = n_double + 1;
                angles(n_double) = theta;
                lengths(n_double) = len;
                for i_mt = 1 : 1 : params.n_mts
                    neighb_found = false;
                    for delta = 1 : 1 : params.max_sites
                        if neighb_found
                            break;
                        end
                        for dir = -1 : 2 : 1
                            i_scan = i_site + (dir * delta);
                            if i_scan >= 1 && i_scan <= params.mt_lengths(i_mt)
                                sid_neighb = occupancy(i_scan, i_mt, i_data + i_window);
                                if sid_neighb == sid_xlink
                                    neighb_distance(n_double + (i_mt - 1)) = delta - 1;
                                    neighb_found = true;
                                    break;
                                end
                            end
                        end
                    end
                end
            % if crosslinker is singly bound, add to statistics
            else
                if pos_x < overlap_start || pos_x > overlap_end
                    n_single_out = n_single_out + 1;
                else
                    n_single_in = n_single_in + 1;
                end
            end
        end
    end
    stats_single_out(i_data) = n_single_out;
    stats_single_in(i_data) = n_single_in;
    stats_double(i_data) = n_double;
    end   
    
    % determine overlap start and end coordinates 
    plus_pos_top = filament_pos(:, 1, 2, i_data);
    minus_pos_top = filament_pos(:, 2, 2, i_data);
    plus_pos_bot = filament_pos(:, 1, 1, i_data);
    minus_pos_bot = filament_pos(:, 2, 1, i_data);
    overlap_start = max(plus_pos_bot(1), minus_pos_top(1));
    overlap_end = min(plus_pos_top(1), minus_pos_bot(1));
    
    
   % disp(max(lengths))
    
    % Plot occupancy profile
    subplot(2, 4, 1)
    occupancy(occupancy ~= sid_xlink) = 0;
    occupancy(occupancy == sid_xlink) = 1;
    avg_top = zeros([params.mt_lengths(2) 1]);
    avg_bot = zeros([params.mt_lengths(1) 1]);   
    for i_window = 0 : 1 : frames_per_plot
        for i_site = 1 : params.mt_lengths(1)
            avg_top(i_site) = avg_top(i_site) + double(occupancy(i_site, 2, i_data + i_window))./ frames_per_plot;
        end
        for i_site = 1 : params.mt_lengths(2)
            avg_bot(i_site) = avg_bot(i_site) + double(occupancy(i_site, 1, i_data + i_window)) ./ frames_per_plot;
        end
    end
    x_begin = minus_pos_top(1);
    x_end = minus_pos_bot(1);
    span = x_end - x_begin;
    span_sites = int32(span / (params.site_size * 1000));
   % i_offset = int32((minus_pos_top(1) - plus_pos_bot(1)) / (params.site_size * 1000)) + 1; 
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
    xlim([(overlap_start - 200) (overlap_end + 200)]);
    ylim([(plus_pos_bot(2) - 2.5) (plus_pos_top(2) + 2.5)]);
    line([plus_pos_top(1), minus_pos_top(1)],[plus_pos_top(2), minus_pos_top(2)], ...
        'LineWidth', 4, 'Color', [0.7 0.7 0.7]);
    line([plus_pos_bot(1), minus_pos_bot(1)],[plus_pos_bot(2), minus_pos_bot(2)], ...
        'LineWidth', 4, 'Color', [0.7 0.7 0.7]);
    for i_site = 1 : 1 : params.max_sites
        id = protein_ids(i_site, 1, i_data);
        sid = occupancy(i_site, 1, i_data);
        if sid == sid_xlink
            n_sites = params.mt_lengths(1);
            bot_mt_vec = [minus_pos_bot(1) - plus_pos_bot(1), minus_pos_bot(2) - plus_pos_bot(2)];
            pos_x = plus_pos_bot(1) + ((i_site-1)/(n_sites-1))*bot_mt_vec(1);
            pos_y = plus_pos_bot(2) + ((i_site-1)/(n_sites-1))*bot_mt_vec(2);
            if partner_indices(i_site, 1, i_data) ~= -1           
                ii_site = partner_indices(i_site, 1, i_data);
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
        id = protein_ids(i_site, 1, i_data);
        sid = occupancy(i_site, 1, i_data);
        if sid == sid_xlink
            n_sites = params.mt_lengths(1);
            bot_mt_vec = [minus_pos_bot(1) - plus_pos_bot(1), minus_pos_bot(2) - plus_pos_bot(2)];
            pos_x = plus_pos_bot(1) + ((i_site-1)/(n_sites-1))*bot_mt_vec(1);
            pos_y = plus_pos_bot(2) + ((i_site-1)/(n_sites-1))*bot_mt_vec(2);
            if partner_indices(i_site, 1, i_data) ~= -1           
                ii_site = partner_indices(i_site, 1, i_data);
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
    x_axis = [start_frame:1:end_frame-1];
    hold all
    plot(x_axis, stats_double, 'o', 'LineWidth', 2);
    plot(x_axis, stats_single_in, 'o', 'LineWidth', 2);
    plot(x_axis, stats_single_out, 'o', 'LineWidth', 2);
    legend(["Double", "Single", "Escaped"], 'location', 'northeast');
    ylim([0 inf]);
    xlim([start_frame end_frame]);
    ylabel("Count (global)");
    xlabel("Timestep number");
    
    % Plot histogram of xlink angles
    subplot(2, 4, 5)
    histfit(angles);
    xlim([35 145]);
    xticks([40 90 140]);
    ylabel("Count");
    xlabel("Crosslinker angle (deg)");
    
    % Plot histogram of xlink lengths
    subplot(2, 4, 6)
    histfit(lengths);
    xlim([15 50]);
    xticks([20 32 45]);
    ylabel("Count");
    xlabel("Crosslinker length (nm)");
    
   % Plot histogram of xlink neighbor distance
    subplot(2, 4, 7)
    n_bins = int32(sqrt(length(neighb_distance)));
    h = histfit(neighb_distance, n_bins, 'exponential');
    xlim([-1 15]);
    xticks([0 4 8]);
    ylabel("Count");
    xlabel("Distance to neighbor (n sites)");
    
    % n_bound histogram here 
    subplot(2, 4, 8)
    bar([n_double n_single_in n_single_out]');
    ylabel("Count");
    xlabel("Populations");
    xticklabels(["Double", "Single", "Escaped"]);

    
    time = (i_data - start_frame) * params.time_per_datapoint;
    str = sprintf('Time: %#.2f seconds', time);
    dim = [0.15 0.69 .3 .3];
    annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on'); 
    str = sprintf('Force: %#.2f pN', xlink_data(1, i_data));
    dim = [0.35 0.69 0.3 0.3];
    annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on'); 
    drawnow();
       writeVideo(v, getframe(gcf));
end

close(v);
