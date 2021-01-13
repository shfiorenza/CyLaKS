clear variables;

fileDirectory = '/home/shane/projects/CyLaKS/%s';

sim_name = 'testo';

movie_name = 'test';
start_frame = 1;
frames_per_plot = 10;
movie_duration = 30; % in seconds

sid_site = 0;
sid_xlink = 1;
sid_motor = 2;

% Colors
blue = [30 144 255] / 255;
purple = [128 0 128] / 255;
color = [purple; blue];

% Open log file and parse it into param labels & their values
log_file = sprintf(fileDirectory, sprintf('%s.log', sim_name));
log = textscan(fileread(log_file), '%s %s', 'Delimiter', '=');
params = log{1, 1};
values = log{1, 2};
% Read in number of MTs
n_mts = sscanf(values{contains(params, "count ")}, '%g');
mt_lengths = zeros(1, n_mts);
for i_mt = 1 : n_mts
    string = sprintf("n_sites[%i] ", i_mt - 1);
    mt_lengths(i_mt) = sscanf(values{contains(params, string)}, '%i');
    if any(contains(params, sprintf("N_SITES[%i] ", i_mt - 1)) ~= 0)
        string = sprintf("N_SITES[%i] ", i_mt - 1);
        mt_lengths(i_mt) = sscanf(values{contains(params, string)}, '%i');
    end
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
site_size = 0.0082; % in um
max_sites = max(mt_lengths);
xlink_cutoff = 5;
teth_cutoff = 19;

r_prot = (site_size*1000);

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
occupancyFileName = '%s_occupancy.file';
partnerFileName = '%s_partner_index.file'; % not needed?
tethFileName = '%s_tether_anchor_pos.file'; % not needed?
motorHeadFileName = '%s_motor_head_trailing.file';
filamentFile = sprintf(fileDirectory, sprintf(filamentFileName, sim_name));
proteinFile = sprintf(fileDirectory, sprintf(proteinFileName, sim_name));
occupancyFile = sprintf(fileDirectory, sprintf(occupancyFileName, sim_name));
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
occupancy = zeros(max_sites, n_mts, n_datapoints) - 1;
if isfile(occupancyFile)
    file = fopen(occupancyFile);
    data = fread(file, n_mts * max_sites * n_datapoints, '*int');
    fclose(file);
    occupancy = reshape(data, max_sites, n_mts, n_datapoints);
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
    file = fopen(motorHeadFile);
    data = fread(file, [n_mts * max_sites * n_datapoints], '*bool');
    fclose(file);
    motor_trailing = reshape(data, max_sites, n_mts, n_datapoints);
    
end

% Run through all datapoints; each one is a frame in our movie
for i_data = start_frame : frames_per_plot : end_frame
    % Clear figure so that it only displays figures from current datapoint
    clf;
    % Set Axes properties
    ax = axes('Units', 'normalized', 'Position', [0.075 0.085 0.9 0.9]);
    hold all;
    min_x = min(min(filament_pos(1, :, :, i_data)));
    max_x = max(max(filament_pos(1, :, :, i_data)));
    min_y = min(min(filament_pos(2, :, :, i_data)));
    max_y = max(max(filament_pos(2, :, :, i_data)));
    y_avg = (min_y + max_y)/2;
    width = (max_x - min_x) + 50;
    height = (3/5 * width); 
    ax.XLim = [(min_x - 25) (max_x + 25)];
    ax.YLim = [y_avg - height/2 y_avg + height/2];
    %ax.XLim = [-167 167];
    %ax.YLim = [-100 100];
    %ax.XTick = linspace(roundn(min_x, 2), roundn(max_x, 2), 5);
    %ax.YTick = linspace(roundn((3/5)*min_x, 2), roundn((3/5)*max_x, 2), 5);
    ax.TickLength = [0.005 0.005];
    ax.XLabel.String = 'x position (nm)';
    ax.YLabel.String = 'y position (nm)';
    % Draw filaments
    if(n_mts > 1)
        com_y_one = (filament_pos(2, 1, 1, i_data) + filament_pos(2, 2, 1, i_data))/2;
        com_y_two = (filament_pos(2, 1, 2, i_data) + filament_pos(2, 2, 2, i_data))/2;
        disp(com_y_one - com_y_two);
    end
    for i_mt = 1:1:n_mts
        plus_pos = filament_pos(:, 1, i_mt, i_data);
        minus_pos = filament_pos(:, 2, i_mt, i_data);
        line([plus_pos(1)-r_prot/2, minus_pos(1)-r_prot/2],[plus_pos(2), minus_pos(2)], ...
            'LineWidth', height / 100, 'Color', [0.7 0.7 0.7]);
        n_sites = mt_lengths(i_mt);
        % Draw proteins
        for i_site = 1 : n_sites
            id = protein_ids(i_site, i_mt, i_data);
            sid = occupancy(i_site, i_mt, i_data);
            if(id ~= -1)
                if i_mt == 1
                    dx = -1;
                    mt_dir = 1;
                    line_vec = [minus_pos(1) - plus_pos(1), minus_pos(2) - plus_pos(2)];
                    pos_x = plus_pos(1) + ((i_site-1)/n_sites)*line_vec(1);
                    pos_y = plus_pos(2) + ((i_site-1)/n_sites)*line_vec(2);
                else
                    dx = 1;
                    mt_dir = -1;
                    line_vec = [plus_pos(1) - minus_pos(1), plus_pos(2) - minus_pos(2)];
                    pos_x = minus_pos(1) + ((i_site-1)/n_sites)*line_vec(1);
                    pos_y = minus_pos(2) + ((i_site-1)/n_sites)*line_vec(2);
                end
                if sid == sid_xlink
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
                            line([pos_x, endpos_x],[pos_y, endpos_y], ...
                                'LineWidth', 1, 'Color', purple);
                        end
                    end
                elseif sid == sid_motor
                    singly_bound = true;
                    i_fwd = i_site + 1;
                    if i_fwd <= n_sites
                        if protein_ids(i_fwd, i_mt, i_data) == id
                            endpos_x = pos_x + r_prot/2;
                            endpos_y = pos_y + mt_dir * r_prot;
                            line([pos_x, endpos_x], [pos_y, endpos_y], ...
                                'LineWidth', 1, 'Color', blue);
                            singly_bound = false;
                        end
                    end
                    i_bck = i_site - 1;
                    if i_bck > 0
                        if protein_ids(i_bck, i_mt, i_data) == id
                            endpos_x = pos_x - r_prot/2;
                            endpos_y = pos_y + mt_dir  * r_prot;
                            line([pos_x, endpos_x], [pos_y, endpos_y], ...
                                'LineWidth', 1, 'Color', blue);
                            singly_bound = false;
                        end
                    end
                    if singly_bound
                        endpos_x = pos_x;
                        endpos_y = pos_y + mt_dir * 1.12 * r_prot;
                        line([pos_x, endpos_x], [pos_y, endpos_y], ...
                            'LineWidth', 1, 'Color', blue);
                        dir = -dx; 
                        if motor_trailing(i_site, i_mt, i_data)
                           dir = dx; 
                        end
                        dockpos_x = pos_x + 1.045 * r_prot * dir;
                        dockpos_y = endpos_y - mt_dir * r_prot / 4;
                        line([endpos_x, dockpos_x], [endpos_y, dockpos_y], ...
                            'LineWidth', 1, 'Color', blue);
                        rectangle('Position', [dockpos_x - r_prot/2 dockpos_y - r_prot/2 r_prot r_prot], ...
                            'FaceColor', color(sid, :), 'Curvature', [1 1]);
                    end
                end
                % Draw protein head
                pos = [pos_x-(r_prot/2) pos_y-(r_prot/2) r_prot r_prot];
                rectangle('Position', pos,'FaceColor', color(sid, :), 'Curvature', [1 1]);
            end 
        end
    end
    dim = [0.11 0.625 .3 .3];
    time = (i_data - start_frame) * time_per_datapoint;
    str = sprintf('Time: %#.2f seconds', time);
    annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on');
    drawnow();
    %  frame = getframe(gcf); %(fig1); %, frame_box);
    writeVideo(v, getframe(gcf));
end

close(v);
%}