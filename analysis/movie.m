clear variables;
sim_name = 'test'; % Raw sim name; do not include directory
movie_name = 'testin2';

% Movie details
start_frame = 1;
frames_per_plot = 10;
movie_duration = 30; % in seconds

% Species IDs
sid_site = 0;
sid_xlink = 1;
sid_motor = 2;

% Colors
blue = [30 144 255] / 255;
purple = [128 0 128] / 255;
color = [purple; blue];

% Open log file and parse it into param labels & their values
fileDirectory = '../%s';
%fileDirectory = '/home/shane/data_kif4a_paper/run_mobility_both/%s';
log_file = sprintf(fileDirectory, sprintf('%s.log', sim_name));
log = textscan(fileread(log_file), '%s %s', 'Delimiter', '=');
params = log{1, 1};
values = log{1, 2};
% Read in system params
dt = sscanf(values{contains(params, "dt ")}, '%g');
time_per_datapoint = sscanf(values{contains(params, "t_snapshot ")}, '%g');
n_datapoints = str2double(values{contains(params, "n_datapoints ")});
% Use actual recorded number of datapoints to parse thru data/etc
if any(contains(params, "N_DATAPOINTS ") ~= 0)
    n_datapoints = str2double(values{contains(params, "N_DATAPOINTS ")});
end
site_size =  sscanf(values{contains(params, "site_size ")}, '%g') / 1000; % in um
% Read in number of MTs
n_mts = sscanf(values{contains(params, "count ")}, '%g');
if any(contains(params, "COUNT ") ~= 0)
    n_mts = sscanf(values{contains(params, "COUNT ")}, '%g');
end
if any(contains(params, "n_subfilaments") ~= 0)
    n_sub = sscanf(values{contains(params, "n_subfilaments ")}, '%g');
    if n_sub > n_mts
       n_mts = n_sub;
    end
end
% Read in MT lengths (in n_sites)
mt_lengths = zeros(1, n_mts);
for i_mt = 1 : n_mts
    string = sprintf("n_sites[%i] ", i_mt - 1);
    mt_lengths(i_mt) =  sscanf(values{contains(params, string)}, '%i');
    if any(contains(params, sprintf("N_SITES[%i] ", i_mt - 1)) ~= 0)
        string = sprintf("N_SITES[%i] ", i_mt - 1);
        mt_lengths(i_mt) = sscanf(values{contains(params, string)}, '%i');
    end
end
max_sites = max(mt_lengths);
polarity = zeros(1, n_mts);
for i_mt = 1 : n_mts
    string = sprintf("polarity[%i] ", i_mt - 1);
    polarity(i_mt) =  sscanf(values{contains(params, string)}, '%i');
end


n_dims = 2; % hard-coded for now; CyLaKS always outputs data in 2-D
xlink_cutoff = 5; % FIXME: dynamically get this from log
teth_cutoff = 10; % FIXME: dynamically get this from log

r_prot = (site_size*1000);

end_frame = n_datapoints;
% Figure parameters (i.e., how they appear)
site_height = 1;
site_width = 1;
active_frames = end_frame - start_frame;

% Videowriter details
v = VideoWriter(movie_name);
v.FrameRate = (active_frames / frames_per_plot) / movie_duration;
open(v);
frame_box = [0 0 1445 200];

% Figure details
fig1 = figure;
set(fig1, 'Position', [50 50 1000 500]);

% File info
filamentFileName = '%s_filament_pos.file';
proteinFileName = '%s_protein_id.file';
occupancyFileName = '%s_occupancy.file';
motorHeadFileName = '%s_motor_head_trailing.file';

partnerFileName = '%s_partner_index.file'; % not needed?
tethFileName = '%s_tether_anchor_pos.file'; % not needed?

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
    file = fopen(tethFile);
    data = fread(file, [n_mts * max_sites * n_datapoints], '*double');
    fclose(file);
    teth_data = reshape(data, max_sites, n_mts, n_datapoints);
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
    avg_y = (min_y + max_y)/2;
    height = 400; % (1/10)*(max_x - min_x);
    ax.XLim = [(min_x - 25) (max_x + 25)];
    %ax.XLim = [(min_x - 25) (min_x + 475)];
    ax.YLim = [(avg_y - height/2) (avg_y + height/2)];
    ax.XTick = linspace(roundn(min_x, 2), roundn(max_x, 2), 5);
    ax.YTick = linspace(roundn(avg_y - height/2, 2), roundn(avg_y + height/2, 2), 3);
    ax.TickLength = [0.005 0.005];
    ax.XLabel.String = 'x position (nm)';
    ax.YLabel.String = 'y position (nm)';
    % Draw filaments
    if(n_mts > 1)
        com_y_one = (filament_pos(2, 1, 1, i_data) + filament_pos(2, 2, 1, i_data))/2;
        com_y_two = (filament_pos(2, 1, 2, i_data) + filament_pos(2, 2, 2, i_data))/2;
        %disp(com_y_one - com_y_two);
    end
    for i_mt = 1:1:n_mts
        plus_pos = filament_pos(:, 1, i_mt, i_data);
        minus_pos = filament_pos(:, 2, i_mt, i_data);
        line([plus_pos(1)-r_prot/2, minus_pos(1)-r_prot/2],[plus_pos(2), minus_pos(2)], ...
            'LineWidth', 4, 'Color', [0.7 0.7 0.7]);
        n_sites = mt_lengths(i_mt);
        dx = -1;
        mt_dir = 1;
        line_vec = [minus_pos(1) - plus_pos(1), minus_pos(2) - plus_pos(2)];
        if polarity(i_mt) == 1
            dx = 1;
            mt_dir = -1;
            line_vec = [plus_pos(1) - minus_pos(1), plus_pos(2) - minus_pos(2)];
        end
        
        
        % Draw proteins
        for i_site = 1 : n_sites
            id = protein_ids(i_site, i_mt, i_data);
            sid = occupancy(i_site, i_mt, i_data);
            if(id ~= -1)
                pos_x = plus_pos(1) + ((i_site-1)/(n_sites-1))*line_vec(1);
                pos_y = plus_pos(2) + ((i_site-1)/(n_sites-1))*line_vec(2);
                if polarity(i_mt) == 1
                    pos_x = minus_pos(1) + ((i_site-1)/(n_sites-1))*line_vec(1);
                    pos_y = minus_pos(2) + ((i_site-1)/(n_sites-1))*line_vec(2);
                end
                
                if sid == sid_xlink
                    % Draw spring connecting crosslinker if appropriate
                    if(i_mt == 0) % lol
                        if(partner_indices(i_site, i_mt, i_data) ~= -1)
                            ii_site = partner_indices(i_site, i_mt, i_data);
                            nn_sites = mt_lengths(2);
                            p_pos = filament_pos(:, 1, 2, i_data);
                            m_pos = filament_pos(:, 2, 2, i_data);
                            neighb_vec = [p_pos(1) - m_pos(1), p_pos(2) - m_pos(2)];
                            endpos_x = m_pos(1) + (double(ii_site)/(nn_sites-1))*neighb_vec(1);
                            endpos_y = m_pos(2) + (double(ii_site)/(nn_sites-1))*neighb_vec(2);
  %                          line([pos_x, endpos_x],[pos_y, endpos_y], ...
  %                              'LineWidth', 1, 'Color', purple);
                        else
  %                          line([pos_x, pos_x], [pos_y, pos_y + 3* r_prot], ...
  %                              'LineWidth', 1, 'Color', purple);
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
        % Draw tethers
        
        line_vec = [minus_pos(1) - plus_pos(1), minus_pos(2) - plus_pos(2)];
        teth_coords = teth_data(:, i_mt, i_data);
        for i_teth=1:1:max_sites - 1
            if(teth_coords(i_teth) ~= -1)
                start_x = plus_pos(1) + ((i_teth-1)/(n_sites-1))*line_vec(1);
                start_y = plus_pos(2) + ((i_teth-1)/(n_sites-1))*line_vec(2) + 1.045 * r_prot;
                end_x = teth_coords(i_teth); % - 4.1;
                end_y = start_y + 16;
               % xa = start_x; ya = start_y;
               % xb = end_x; yb = end_y + 50;
               % ne = 5; a = 2; ro = 50;
                line([start_x, end_x], [start_y, end_y], ...
                    'LineWidth', 1, 'Color', blue);
                %[xs,ys] = spring(xa,ya,xb,yb,ne,a,ro);
                %plot(xs,ys,'LineWidth', 1, 'Color', 'black');
                
                %disp(teth_coords(i_teth))
                %{
                end_height = (mt_height + neighb_mt_height + site_height) / 2;
                if(n_mts < 2)
                    end_height = mt_height + 5*site_height; 
                end
                start_height = mt_height + 5*site_height / 2;
                if(mod(i_mt, 2) == 0)
                    start_height = mt_height - 3*site_height / 2;
                end
               
                if(teth_coords(i_teth) ~= teth_coords(i_teth + 1))
                    start_pos = i_teth*site_width + mt_pos - 1;
                    if(i_teth > 1)
                        if(motor_IDs(i_teth) ~= motor_IDs(i_teth - 1) ...
                        && motor_IDs(i_teth) ~= motor_IDs(i_teth + 1))
                            start_pos = start_pos + site_width/2 - 1;
                        end
                    elseif(motor_IDs(i_teth) ~= motor_IDs(i_teth + 1))
                        start_pos = start_pos + site_width/2 - 1;
                    end                    
                    end_pos = teth_coords(i_teth)*site_width + (3/2)*site_width - 1;
                    xa = start_pos; ya = start_height;
                    xb = end_pos; yb = end_height;
                    ne = 10; a = 10; ro = 0.5;
                    if abs(xa - xb) <= cutoff
                        [xs,ys] = spring(xa,ya,xb,yb,ne,a,ro);
                        plot(xs,ys,'LineWidth', 1, 'Color', 'black');
                    else
                        disp(xa - xb);
                        disp(teth_coords(i_teth));
                        disp(i_teth);
                    end
               end
                %}
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
