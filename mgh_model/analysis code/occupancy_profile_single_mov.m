clear variables;
% Often-changed variables
off_ratio = 10;
k_hydrolyze = 90;
jam_ratio = 700;
n_sites = 1000;
% Pseudo-constant variables
motor_speciesID = 2;
xlink_speciesID = 1;
n_steps = 1000000;
n_datapoints = 10000;
steps_per_plot = 1000;
starting_point = 1;
active_datapoints = n_datapoints - starting_point;
delta_t = 0.0001;
time_per_frame = delta_t * (n_steps / n_datapoints);

%simName = sprintf('Endtag_%i', n_sites);
if(off_ratio == 1)
    fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/scan_output/%s';
    simName = sprintf('k_hydrolyze_%i/jam_ratio_%i/endtag_%i_%i_%i', ...
        k_hydrolyze, jam_ratio, k_hydrolyze, jam_ratio, n_sites);
else
   fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
   % simName = sprintf('off_ratio_%i/k_hydrolyze_%i/jam_ratio_%i/endtag_%i_%i_%i_%i', ...
    %    off_ratio, k_hydrolyze, jam_ratio, off_ratio, k_hydrolyze, jam_ratio, n_sites);
    simName = sprintf('endtag_%i_%i_%i_%i', off_ratio, k_hydrolyze, jam_ratio, n_sites);
end
fileStruct = '%s_occupancy.file';
legendLabel = {'Motors', 'Crosslinkers', 'Combined'};

%movie_name = sprintf('%s_mov.avi', simName);
movie_name = 'test.avi';
% Videowriter details
v = VideoWriter(movie_name);
v.FrameRate = (active_datapoints / steps_per_plot) / 30;
open(v);
frame_box = [0, 0, 1.5*480, 1.5*300];

fileName = sprintf(fileDirectory, sprintf(fileStruct, simName));
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
        motor_smoothed_avg = smoothdata(motor_avg_occupancy);
        xlink_smoothed_avg = smoothdata(xlink_avg_occupancy);
        net_smoothed_avg = motor_smoothed_avg + xlink_smoothed_avg;
        occupancy_slope = abs(gradient(net_smoothed_avg, 0.08));
        max_slope = max(occupancy_slope);
        
        motor_avg_occupancy = zeros([n_sites 1]);
        xlink_avg_occupancy = zeros([n_sites 1]);
        
        endtag_site = 0;
        past_max = false;
        for i_site=1:n_sites
            if(occupancy_slope(i_site) >= max_slope)
                past_max = true;
            end
            if(occupancy_slope(i_site) < max_slope/2 && past_max)
                endtag_site = i_site;
                break;
            end
        end
        endtag_length = endtag_site*0.008;
        
        %%plot fig%%

        clf;
        ax = axes('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
        hold all

        plot(linspace(0, n_sites*0.008, n_sites), motor_smoothed_avg);
        
       % plot(linspace(0, n_sites*0.008, n_sites), occupancy_slope);
        %plot(linspace(0, n_sites*0.008, n_sites), xlink_smoothed_avg);
        %plot(linspace(0, n_sites*0.008, n_sites), net_smoothed_avg);
        
        % Put vertical red line where endtag starting position is
        plot([endtag_length endtag_length], [0 1], ':r', 'LineWidth', 0.1);
        
        %%style stuff%%
        
        
        %title({sprintf('%g micron-long microtubule', n_sites*0.008), ...
        %sprintf('%#.1f nM PRC1 and %#.1f nM kinesin-1', xlink_conc, motor_conc), ...
        title(sprintf('Endtag length: %g microns', endtag_length));
        
        xlabel({'Distance along microtubule relative to plus-end (microns)'});
        ylabel('Fraction of the time occupied');
        %xlim([0 0.5]);
        ylim([0 1]);
        
        %grid on
        %grid minor
        axis = gca;
        axis.TickDir = 'out';
        axis.Box = 'off';
        axis.GridLineStyle = '-';
        set(findall(axis, 'Type', 'Line'), 'LineWidth', 2);
        %legend(legendLabel, 'Location', 'northeast');
       
    
    dim = [0.7 0.57 .3 .3];
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