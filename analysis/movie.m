clear variables;

sim_name = 'shep_multiPF_0_0.131_4'; % Raw sim name; do not include directory
sim_name = 'test10x175';
output_movie_name = 'test';

start_frame = 1;
end_frame = -1;  % set to -1 to run until end of data

frames_per_plot = 10; 
movie_duration = 30; % in seconds

% Load parameter structure
file_dir = '..';  % Default; only change if you move CyLaKS output files
params = load_parameters(sprintf('%s/%s', file_dir, sim_name));

% Set plotting helper variables
sid_site = 0;
sid_xlink = 1;
sid_motor = 2;
r_prot = (params.site_size*1000);
site_height = 1;
site_width = 1;
blue = [30 144 255] / 255;
purple = [128 0 128] / 255;
color = [purple; blue];
if end_frame == -1
    end_frame = params.n_datapoints; 
end
active_frames = end_frame - start_frame;

% Initialize videowriter object
v = VideoWriter(output_movie_name);
v.FrameRate = (active_frames / frames_per_plot) / movie_duration;
open(v);
frame_box = [0 0 1445 200];

% Open figure and set to desired size (each frame must be this same size)
fig1 = figure;
set(fig1, 'Position', [50 50 1000 500]);

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

motor_trailing_filename = sprintf('%s/%s_motor_head_trailing.file', file_dir, sim_name);
motor_trailing = zeros(params.max_sites, params.n_mts, params.n_datapoints);
%motor_trailing = load_data(motor_trailing, motor_trailing_filename, '*bool');

partner_filename = sprintf('%s/%s_partner_index.file', file_dir, sim_name);
partner_indices = zeros(params.max_sites, params.n_mts, params.n_datapoints) - 1;
partner_indices = load_data(partner_indices, partner_filename, '*int');

teth_filename = sprintf('%s/%s_tether_anchor_pos.file', file_dir, sim_name);
teth_data = zeros(params.max_sites, params.n_mts, params.n_datapoints) - 1;
teth_data = load_data(teth_data, teth_filename, '*double');

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
    %top_start = filament_pos(1, 2, 2, i_data);
    %top_end = filament_pos(1, 1, 2, i_data);
    ax.XLim = [(min_x - 25) (max_x + 25)];
    %ax.XLim = [(min_x - 25) (min_x + 475)];
    %ax.XLim = [(top_start - 25) (top_end + 25)];
    ax.YLim = [(avg_y - height/2) (avg_y + height/2)];
    %ax.XTick = linspace(roundn(min_x, 2), roundn(max_x, 2), 5);
    ax.XTick = linspace(-2000, 2000, 11);
    ax.YTick = linspace(roundn(avg_y - height/2, 2), roundn(avg_y + height/2, 2), 3);
    ax.TickLength = [0.005 0.005];
    ax.XLabel.String = 'x position (nm)';
    ax.YLabel.String = 'y position (nm)';
    % Draw filaments
    if(params.n_mts > 1)
        com_y_one = (filament_pos(2, 1, 1, i_data) + filament_pos(2, 2, 1, i_data))/2;
        com_y_two = (filament_pos(2, 1, 2, i_data) + filament_pos(2, 2, 2, i_data))/2;
    end
    for i_mt = 1:1:params.n_mts
        plus_pos = filament_pos(:, 1, i_mt, i_data);
        minus_pos = filament_pos(:, 2, i_mt, i_data);
        line([plus_pos(1)-r_prot/2, minus_pos(1)-r_prot/2],[plus_pos(2), minus_pos(2)], ...
            'LineWidth', 4, 'Color', [0.7 0.7 0.7]);
        n_sites = params.mt_lengths(i_mt);
        dx = -1;
        mt_dir = 1;
        line_vec = [minus_pos(1) - plus_pos(1), minus_pos(2) - plus_pos(2)];
        if params.polarity(i_mt) == 1
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
                if params.polarity(i_mt) == 1
                    pos_x = minus_pos(1) + ((i_site-1)/(n_sites-1))*line_vec(1);
                    pos_y = minus_pos(2) + ((i_site-1)/(n_sites-1))*line_vec(2);
                end
                
                if sid == sid_xlink
                    % Draw spring connecting crosslinker if appropriate
                    if(i_mt == 1) % lol
                        if(partner_indices(i_site, i_mt, i_data) ~= -1)
                            ii_site = partner_indices(i_site, i_mt, i_data);
                            nn_sites = params.mt_lengths(2);
                            p_pos = filament_pos(:, 1, 2, i_data);
                            m_pos = filament_pos(:, 2, 2, i_data);
                            neighb_vec = [p_pos(1) - m_pos(1), p_pos(2) - m_pos(2)];
                            endpos_x = m_pos(1) + (double(ii_site)/(nn_sites-1))*neighb_vec(1);
                            endpos_y = m_pos(2) + (double(ii_site)/(nn_sites-1))*neighb_vec(2);
                            line([pos_x, endpos_x],[pos_y, endpos_y], ...
                                'LineWidth', 1, 'Color', purple);
                        else
                           % line([pos_x, pos_x], [pos_y, pos_y + 3* r_prot], ...
                           %     'LineWidth', 1, 'Color', purple);
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
        for i_teth=1:1:params.max_sites - 1
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
                if(params.n_mts < 2)
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
    time = (i_data - start_frame) * params.time_per_datapoint;
    str = sprintf('Time: %#.2f seconds', time);
    annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on');
    drawnow();
    %  frame = getframe(gcf); %(fig1); %, frame_box);
    writeVideo(v, getframe(gcf));
end

close(v);
%}
