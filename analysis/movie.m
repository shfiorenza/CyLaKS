clear variables;

fileDirectory = '/home/shane/projects/CyLaKS/%s';

sim_name = 'test';

movie_name = 'test';
start_frame = 1;
frames_per_plot = 10;
movie_duration = 30; % in seconds

r_prot = 10;

% Open log file and parse it into param labels & their values
log_file = sprintf(fileDirectory, sprintf('%s.log', sim_name));
log = textscan(fileread(log_file), '%s %s', 'Delimiter', '=');
params = log{1, 1};
values = log{1, 2};
% Read in number of MTs
n_mts = str2double(values{contains(params, "count ")});
mt_lengths = zeros(1, n_mts);
for i_mt = 1 : n_mts
    string = sprintf("n_sites[%i] ", i_mt - 1);
    mt_lengths(i_mt) = sscanf(values{contains(params, string)}, '%i');
end
% Read in system params
dt = sscanf(values{contains(params, "dt ")}, '%g');
steps_per_datapoint = str2double(values{contains(params, "n_steps_per_snapshot ")});
time_per_datapoint = dt * steps_per_datapoint;
n_datapoints = str2double(values{contains(params, "n_datapoints ")});
% Use actual recorded number of datapoints to parse thru data/etc
if any(contains(params, "N_DATAPOINTS ") ~= 0)
    n_datapoints = str2double(values{contains(params, "N_DATAPOINTS ")});
end
n_dims = 2;
site_size = 0.008; % in um
max_sites = max(mt_lengths);
xlink_cutoff = 5;
teth_cutoff = 19;

% Colors
blue = [30 144 255] / 255;
purple = [128 0 128] / 255;

% Figure parameters (i.e., how they appear)
site_height = 1;
site_width = 1;
end_frame = n_datapoints;
active_frames = end_frame - start_frame;

% Videowriter details
v = VideoWriter(movie_name);
v.FrameRate = (active_frames / frames_per_plot) / movie_duration;
open(v);
frame_box = [0 0 1445 200];

% Figure details
fig1 = figure;
set(fig1, 'Position', [50 50 1000 600]);

% File info
filamentFileName = '%s_filament_pos.file';
proteinFileName = '%s_protein_id.file';
partnerFileName = '%s_partner_index.file';
tethFileName = '%s_tether_anchor_pos.file';
motorHeadFileName = '%s_motor_head_trailing.file';
filamentFile = sprintf(fileDirectory, sprintf(filamentFileName, sim_name));
proteinFile = sprintf(fileDirectory, sprintf(proteinFileName, sim_name));
partnerFile = sprintf(fileDirectory, sprintf(partnerFileName, sim_name));
tethFile = sprintf(fileDirectory, sprintf(tethFileName, sim_name));
motorHeadFile = sprintf(fileDirectory, sprintf(motorHeadFileName, sim_name));

filament_pos = zeros(n_dims, 2, n_mts, n_datapoints);
if isfile(filamentFile)
    file = fopen(filamentFile);
    data = fread(file, 2*n_dims * n_mts * n_datapoints, '*double');
    fclose(file);
    filament_pos = reshape(data, n_dims, 2, n_mts, n_datapoints);
end
protein_ids = zeros(max_sites, n_mts, n_datapoints) - 1;
if isfile(proteinFile)
    file = fopen(proteinFile);
    data = fread(file, n_mts * max_sites * n_datapoints, '*int');
    fclose(file);
    protein_ids = reshape(data, max_sites, n_mts, n_datapoints);
end
partner_indices = zeros(max_sites, n_mts, n_datapoints) - 1;
if isfile(partnerFile)
    file = fopen(partnerFile);
    data = fread(file, n_mts * max_sites * n_datapoints, '*int');
    fclose(file);
    partner_indices = reshape(data, max_sites, n_mts, n_datapoints);
end

teth_data = zeros(max_sites, n_mts, n_datapoints) - 1;
if isfile(tethFile)
    %{
    file = fopen(tethFile);
    data = fread(file, [n_mts * max_sites * n_datapoints], '*double');
    fclose(file);
    teth_data = reshape(data, max_sites, n_mts, n_datapoints);
    %}
end
motor_trailing = zeros(max_sites, n_mts, n_datapoints);
if isfile(motorHeadFile)
    %{
    file = fopen(motorHeadFile);
    data = fread(file, [n_mts * max_sites * n_datapoints], '*bool');
    fclose(file);
    motor_trailing = reshape(data, max_sites, n_mts, n_datapoints);
    %}
end

% Run through all datapoints; each one is a frame in our movie
for i_data = start_frame : frames_per_plot : end_frame
    % Clear figure so that it only displays figures from current datapoint
    clf;
    % Set Axes properties
    ax = axes('Units', 'normalized', 'Position', [0.05 0.075 0.9 0.9]);
    hold all;
    min_x = min(min(filament_pos(1, :, :, i_data)));
    max_x = max(max(filament_pos(1, :, :, i_data)));
    ax.XLim = [(min_x - 25) (max_x + 25)];
    ax.YLim = [(3/5)*(min_x - 25) (3/5)*(max_x + 25)];
    ax.XTick = linspace(roundn(min_x, 2), roundn(max_x, 2), 5);
    ax.YTick = linspace(roundn((3/5)*min_x, 2), roundn((3/5)*max_x, 2), 5);
    ax.TickLength = [0.005 0.005];
    ax.XLabel.String = 'Position in vertical dimension (nm)';
    ax.XLabel.String = 'Position in horizontal dimension (nm)';
    % Draw filaments
    for i_mt = 1:1:n_mts
        plus_pos = filament_pos(:, 1, i_mt, i_data);
        minus_pos = filament_pos(:, 2, i_mt, i_data);
        line([plus_pos(1), minus_pos(1)],[plus_pos(2), minus_pos(2)], 'LineWidth', 4);
        %length = sqrt((plus_pos(1) - minus_pos(1))^2 + (plus_pos(2) - minus_pos(2))^2)
        n_sites = mt_lengths(i_mt);
        % Draw proteins
        for i_site = 1 : n_sites
            if(protein_ids(i_site, i_mt, i_data) ~= -1)
                if i_mt == 1
                    line_vec = [minus_pos(1) - plus_pos(1), minus_pos(2) - plus_pos(2)];
                    pos_x = plus_pos(1) + ((i_site-1)/n_sites)*line_vec(1);
                    pos_y = plus_pos(2) + ((i_site-1)/n_sites)*line_vec(2);
                else
                    line_vec = [plus_pos(1) - minus_pos(1), plus_pos(2) - minus_pos(2)];
                    pos_x = minus_pos(1) + ((i_site-1)/n_sites)*line_vec(1);
                    pos_y = minus_pos(2) + ((i_site-1)/n_sites)*line_vec(2);
                end
                % Draw spring connecting crosslinker if appropriate
                if(n_mts > 1 && i_mt == 1)
                    if(partner_indices(i_site, i_mt, i_data) ~= -1)
                        ii_site = partner_indices(i_site, i_mt, i_data);
                        nn_sites = mt_lengths(2); 
                        p_pos = filament_pos(:, 1, 2, i_data);
                        m_pos = filament_pos(:, 2, 2, i_data); 
                        neighb_vec = [p_pos(1) - m_pos(1), p_pos(2) - m_pos(2)];
                        endpos_x = m_pos(1) + (double(ii_site)/nn_sites)*neighb_vec(1);
                        endpos_y = m_pos(2) + (double(ii_site)/nn_sites)*neighb_vec(2);
                        %line([pos_x, endpos_x],[pos_y, endpos_y], ...
                        %    'LineWidth', 1, 'Color', purple);
                        ne = 4; a = 2; ro = 4;
                        [xs, ys] = spring(pos_x, pos_y, endpos_x, endpos_y, ne, a, ro);
                        plot(xs, ys, 'LineWidth', 1, 'Color', purple);
                    end
                end
                % Draw protein head
                pos = [pos_x-(r_prot/2) pos_y-(r_prot/2) r_prot r_prot];
                rectangle('Position', pos,'FaceColor', purple, 'Curvature', [1 1]);
            end 
        end
        
        
        %{
        if (n_mts > 1)
            
            if (mod(i_mt, 2) == 0)
                neighb_IDs = xlink_data(:, i_mt - 1, i_data);
                neighb_mt_pos = filament_pos(i_mt - 1, i_data) * site_width;
                neighb_mt_height = 8 * (i_mt - 2) * site_height;
            else
                neighb_IDs = xlink_data(:, i_mt + 1, i_data);
                neighb_mt_pos = filament_pos(i_mt + 1, i_data) * site_width;
                neighb_mt_height = 8 * (i_mt) * site_height;
            end
            
        end
        
        % Scan thru MTs and plot all xlinks
        for i_xlink = 1:1:n_sites
            xlink_pos = mt_pos + i_xlink * site_width - 1;
            xlink_height = mt_height + site_height;
            xlink_center_x = xlink_pos + site_width / 2;
            xlink_center_y = xlink_height + site_height / 2;
            
            if (mod(i_mt, 2) == 0)
                xlink_height = mt_height - site_height;
                xlink_center_y = xlink_height - site_height / 2;
            end
            
            if (xlink_IDs(i_xlink) ~= -1)
                double_bound = false;
                
                for i_scan = -xlink_cutoff:1:xlink_cutoff
                    i_neighb = i_xlink + i_scan + mt_pos - neighb_mt_pos;
                    
                    if i_neighb < 1
                        i_neighb = 1;
                    elseif i_neighb > n_sites
                        i_neighb = n_sites;
                    end
                    
                    if (xlink_IDs(i_xlink) == neighb_IDs(i_neighb))
                        neighb_center_x = neighb_mt_pos + (i_neighb - 1/2) * site_width;
                        neighb_center_y = neighb_mt_height + 3 * site_height / 2;
                        neighb_height = neighb_mt_height + site_height;
                        
                        if (mod(i_mt, 2) ~= 0)
                            neighb_center_y = neighb_mt_height + site_height / 2;
                            neighb_height = neighb_mt_height;
                            % To ensure spring is only plotted once per
                            % xlink, only plot on odd-numbered MTs
                            xa = xlink_center_x; ya = xlink_height + site_height;
                            xb = double(neighb_center_x); yb = neighb_height - site_height;
                            ne = 6; a = 6; ro = 1;
                            [xs, ys] = spring(xa, ya, xb, yb, ne, a, ro);
                            plot(xs, ys, 'LineWidth', 1, 'Color', purple);
                        end
                        
                        double_bound = true;
                        rectangle('Position', [xlink_pos xlink_height site_width site_height], ...
                            'FaceColor', purple, 'Curvature', [0.5 0.5]);
                    end
                    
                end
                
                if (double_bound == false)
                    rectangle('Position', [xlink_pos xlink_height site_width site_height], ...
                        'FaceColor', 'm', 'Curvature', [0.5 0.5]);
                    xa = xlink_center_x; ya = xlink_height + site_height;
                    xb = xlink_center_x; yb = xlink_height + 4 * site_height;
                    
                    if (mod(i_mt, 2) == 0)
                        ya = xlink_height;
                        yb = xlink_height - 3 * site_height;
                    end
                    

                end
                
            end
            
        end
        
        % Array of tether coords for this MT
        teth_coords = teth_data(:, i_mt, i_data);
        
        for i_teth = 1:1:max_sites - 1
            
            if (teth_coords(i_teth) ~= -1)
                end_height = (mt_height + neighb_mt_height + site_height) / 2;
                
                if (n_mts < 2)
                    end_height = mt_height + 5 * site_height;
                end
                
                start_height = mt_height + 5 * site_height / 2;
                
                if (mod(i_mt, 2) == 0)
                    start_height = mt_height - 3 * site_height / 2;
                end
                
                if (teth_coords(i_teth) ~= teth_coords(i_teth + 1))
                    start_pos = i_teth * site_width + mt_pos - 1;
                    
                    if (i_teth > 1)
                        
                        if (motor_IDs(i_teth) ~= motor_IDs(i_teth - 1) ...
                                && motor_IDs(i_teth) ~= motor_IDs(i_teth + 1))
                            start_pos = start_pos + site_width / 2 - 1;
                        end
                        
                    elseif (motor_IDs(i_teth) ~= motor_IDs(i_teth + 1))
                        start_pos = start_pos + site_width / 2 - 1;
                    end
                    
                    end_pos = teth_coords(i_teth) * site_width + (3/2) * site_width - 1;
                    xa = start_pos; ya = start_height;
                    xb = end_pos; yb = end_height;
                    ne = 10; a = 10; ro = 0.5;
                    if abs(xa - xb) <= teth_cutoff
                        [xs, ys] = spring(xa, ya, xb, yb, ne, a, ro);
                        plot(xs, ys, 'LineWidth', 1, 'Color', 'black');
                    else
                        disp(xa - xb);
                        disp(teth_coords(i_teth));
                        disp(i_teth);
                    end
                end
            end
        end
        %}
    end
    dim = [0.0105 0.62 .3 .3];
    time = (i_data - 1) * time_per_datapoint;
    %time = time - 500;
    str = sprintf('Time: %#.2f seconds', time);
    annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on');
    drawnow();
    %  frame = getframe(gcf); %(fig1); %, frame_box);
    writeVideo(v, getframe(gcf));
end

close(v);
%}