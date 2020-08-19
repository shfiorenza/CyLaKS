clear variables;
% Often-changed variables

sim_name = 'run_endtag_short/endtag_short_8_1750_0';
file_dir = '/home/shane/Projects/overlap_analysis/mgh_model';

% Open log file and parse it into param labels & their values
log_file = sprintf('%s/%s.log', file_dir, sim_name);
log = textscan(fileread(log_file),'%s %s', 'Delimiter', '=');
params = log{1,1};
values = log{1,2};
% Read in system params
n_mts = str2double(values{contains(params, "count")});
n_sites = values{contains(params, "length")};
n_sites = sscanf(n_sites, '%i');
delta_t = sscanf(values{contains(params, "delta_t")}, '%g');
total_steps = str2double(values{contains(params, "n_steps")});
data_threshold = sscanf(values{contains(params, "data_threshold")}, '%g');
if any(contains(params, "DATA_THRESHOLD") ~= 0)
   data_threshold = str2double(values{contains(params, "DATA_THRESHOLD")});
end
n_steps = total_steps - data_threshold;
% Use max possible number of datapoints to calculate time_per_datapoint (as is done in Sim)
n_datapoints = str2double(values{contains(params, "n_datapoints")});
time_per_datapoint = delta_t * n_steps / n_datapoints;
% Use actual recorded number of datapoints to parse thru data/etc
if any(contains(params, "N_DATAPOINTS") ~= 0)
   n_datapoints = str2double(values{contains(params, "N_DATAPOINTS")});
end

site_size = 0.008;
% Pseudo-constant variables
motor_speciesID = 2;
xlink_speciesID = 1;
steps_per_plot = 1000;
starting_point = 1;
active_datapoints = n_datapoints - starting_point;
time_per_frame = delta_t * (n_steps / n_datapoints);

fileStruct = '%s_occupancy.file';
legendLabel = {'Motors', 'Crosslinkers', 'Combined'};

movie_name = 'test.avi';
% Videowriter details
v = VideoWriter(movie_name);
v.FrameRate = (active_datapoints / steps_per_plot) / 30;
open(v);
frame_box = [0, 0, 1.5*480, 1.5*300];

fileName = sprintf("%s/%s", file_dir, sprintf(fileStruct, sim_name));
data_file = fopen(fileName);
motor_raw_data = fread(data_file, [n_sites, n_datapoints], '*int');
xlink_raw_data = motor_raw_data;
fclose(data_file);

motor_avg_occupancy = zeros([n_sites 1]);
xlink_avg_occupancy = zeros([n_sites 1]);

motor_raw_data(motor_raw_data ~= motor_speciesID) = 0;
motor_raw_data(motor_raw_data == motor_speciesID) = 1;
xlink_raw_data(xlink_raw_data ~= xlink_speciesID) = 0;
xlink_raw_data(xlink_raw_data == xlink_speciesID) = 1;

fig1 = figure(1);
set(fig1,'Position', [100, 100, 2*480, 2*300])

% Read in and average occupancy data over all datapoints
for i=starting_point:1:n_datapoints
    
    motor_avg_occupancy(:,1) = motor_avg_occupancy(:,1) + double(motor_raw_data(:,i))./steps_per_plot;
    xlink_avg_occupancy(:,1) = xlink_avg_occupancy(:,1) + double(xlink_raw_data(:,i))./steps_per_plot;
    
    if(mod(i, steps_per_plot) == 0)
        
        smooth_window = 32; % should be equivalent to diffraction limit
        motor_occupancy = smoothdata(motor_avg_occupancy, 'movmean', smooth_window);
        xlink_occupancy = smoothdata(xlink_avg_occupancy, 'movmean', smooth_window);
        net_occupancy = motor_occupancy + xlink_occupancy;
        occupancy_slope = smoothdata(gradient(net_occupancy, 0.008),'movmean', smooth_window);
        occupancy_accel = smoothdata(gradient(occupancy_slope, 0.008), 'movmean', smooth_window);
        max_occupancy = max(net_occupancy);
        min_slope = min(occupancy_slope);
        
        motor_avg_occupancy(:,1) = 0;
        xlink_avg_occupancy(:,1) = 0;
        
        past_threshold = false;
        i_threshold = 0;
        endtag_site = 0;
        for i_site=1:n_sites
            if(~past_threshold && net_occupancy(i_site) < 0.5 * max_occupancy)
                past_threshold = true;
                i_threshold = i_site;
            end
            if(past_threshold && occupancy_slope(i_site) > (0.5 * min_slope))
                endtag_site = i_site;
                break;
            end
        end
        if endtag_site == 0
            [max_accel, i_peak] = max(occupancy_accel(i_threshold + 1:n_sites));
            endtag_site = i_threshold + i_peak;
        end
        endtag_length = endtag_site * site_size;
     
        %%plot fig%%
        
        clf;
        ax = axes('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
        hold all
        
        plot(linspace(0, n_sites*0.008, n_sites), motor_occupancy);
        
        plot(linspace(0, n_sites*0.008, n_sites), xlink_occupancy);
        plot(linspace(0, n_sites*0.008, n_sites), net_occupancy);
        % plot(linspace(0, n_sites*0.008, n_sites), occupancy_slope);
        % plot(linspace(0, n_sites*0.008, n_sites), occupancy_accel);
        
        % Put vertical red line where endtag starting position is
        plot([endtag_length endtag_length], [0 1], ':', 'LineWidth', 0.1, 'Color', 'red');
        %plot(xlim, [0 0], ':r', 'LineWidth', 0.1);
        
        %%style stuff%%
        
        
        %title({sprintf('%g micron-long microtubule', n_sites*0.008), ...
        %sprintf('%#.1f nM PRC1 and %#.1f nM kinesin-1', xlink_conc, motor_conc), ...
        title(sprintf('Endtag length: %g microns for 1 nM PRC1', endtag_length));
        
        xlabel({'Distance along microtubule relative to plus-end (microns)'});
        ylabel('Fraction of the time occupied');
        %xlim([0 1]);
        ylim([0 1]);
        
        %grid on
        %grid minor
        axis = gca;
        axis.TickDir = 'out';
        axis.Box = 'off';
        axis.GridLineStyle = '-';
        set(findall(axis, 'Type', 'Line'), 'LineWidth', 2);
        legend(legendLabel, 'Location', 'northeast');
        
        
        dim = [0.7 0.47 .3 .3];
        time = i * time_per_frame;
        %time = time - 500;
        str = sprintf('Time: %#.2f seconds', time);
        annotation('textbox',dim,'String',str,'FitBoxToText','on');
        
        drawnow();
        frame = getframe(gcf);
        writeVideo(v, getframe(gcf));
    end
end
close(v);