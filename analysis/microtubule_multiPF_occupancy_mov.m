clear variables;

sim_name = 'test2';
output_movie_name = 'test_out';

start_frame = 1;
frames_per_plot = 100; % in n_datapoints; number of timesteps per output plot
end_frame = -1;  % set to -1 to run until end of data
movie_duration = 30; % in seconds

% Load parameter structure
file_dir = '..';  % Default; only change if you move CyLaKS output files
params = load_parameters(sprintf('%s/%s', file_dir, sim_name));

% Open occupancy data file 
occupancy_filename = sprintf('%s/%s_occupancy.file', file_dir, sim_name);
occupancy = zeros(params.max_sites, params.n_mts, params.n_datapoints) - 1;
occupancy = load_data(occupancy, occupancy_filename, '*int');

% Videowriter details
v = VideoWriter(output_movie_name);
v.FrameRate = (params.n_datapoints / frames_per_plot) / 15;
open(v);
frame_box = [0, 0, 1.5 * 480, 1.5 * 300];

xlink_speciesID = 1;
motor_speciesID = 2;

xlink_raw_data = occupancy; 
motor_raw_data = occupancy; 

xlink_raw_data(xlink_raw_data ~= xlink_speciesID) = 0;
xlink_raw_data(xlink_raw_data == xlink_speciesID) = 1;
motor_raw_data(motor_raw_data ~= motor_speciesID) = 0;
motor_raw_data(motor_raw_data == motor_speciesID) = 1;

xlink_avg_occupancy = zeros([params.max_sites params.n_mts]);
motor_avg_occupancy = zeros([params.max_sites params.n_mts]);

motor_avg_occupancy_tot = zeros([params.max_sites 1]);
xlink_avg_occupancy_tot = zeros([params.max_sites 1]);

fig1 = figure(1);
set(fig1, 'Position', [50, 50, 1200, 600])

% Read in and average occupancy data over all datapoints
for i = 1:1:int32(params.n_datapoints)
    for i_pf = 1 : 1 : params.n_mts
        motor_avg_occupancy(:, i_pf) = motor_avg_occupancy(:, i_pf) + double(motor_raw_data(:, i_pf, i)) ./ frames_per_plot;
        xlink_avg_occupancy(:, i_pf) = xlink_avg_occupancy(:, i_pf) + double(xlink_raw_data(:, i_pf, i)) ./ frames_per_plot;
        motor_avg_occupancy_tot(:, 1) = motor_avg_occupancy_tot(:, 1) + double(motor_raw_data(:, i_pf, i)) ./ (frames_per_plot * params.n_mts);
        xlink_avg_occupancy_tot(:, 1) = xlink_avg_occupancy_tot(:, 1) + double(xlink_raw_data(:, i_pf, i)) ./ (frames_per_plot * params.n_mts);
    end
    if (mod(i, frames_per_plot) == 0)
        smooth_window = 32; % should be equivalent to diffraction limit
        motor_occupancy = smoothdata(motor_avg_occupancy, 'movmean', smooth_window);
        xlink_occupancy = smoothdata(xlink_avg_occupancy, 'movmean', smooth_window);
        motor_occupancy_tot = smoothdata(motor_avg_occupancy_tot, 'movmean', smooth_window);
        xlink_occupancy_tot = smoothdata(xlink_avg_occupancy_tot, 'movmean', smooth_window);  
        
        % Reset arrays to zero before we start counting again 
        for i_pf = 1 : 1 : params.n_mts
            motor_avg_occupancy(:, i_pf) = 0;
            xlink_avg_occupancy(:, i_pf) = 0;
        end
        motor_avg_occupancy_tot(:, 1) = 0;
        xlink_avg_occupancy_tot(:, 1) = 0;
        
        % GET ENDTAG LENGTH HERE

        %%plot fig%%
        clf;
        ax = axes('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
        hold all

        %plot(linspace(0, params.max_sites * params.site_size, params.max_sites), motor_occupancy);
        plot(linspace(0, params.max_sites * params.site_size, params.max_sites), xlink_occupancy, 'LineWidth', 1.5);
        plot(linspace(0, params.max_sites * params.site_size, params.max_sites), xlink_occupancy_tot, 'LineWidth', 3);
        %plot(linspace(0, n_sites*0.008, n_sites), net_occupancy);
        % plot(linspace(0, n_sites*0.008, n_sites), occupancy_slope);
        % plot(linspace(0, n_sites*0.008, n_sites), occupancy_accel);
        % Put vertical red line where endtag starting position is
        %plot([endtag_length endtag_length], [0 1], ':', 'LineWidth', 2, 'Color', 'red');

        %title(sprintf('End-tag length: %g microns', endtag_length), 'FontSize', 14)
        xlabel('Distance from plus-end (\mum)', 'FontSize', 14);
        ylabel('Fractional site occupancy', 'FontSize', 14);
        ylim([0 1]);
        five_percent = params.max_sites * params.site_size / 20.0;
        xlim([-five_percent params.max_sites * params.site_size + five_percent]);
        set(gca, 'FontSize', 14);
        axis = gca;
        %axis.TickDir = 'out';
        axis.Box = 'off';
        axis.GridLineStyle = '-';
        %set(findall(axis, 'Type', 'Line'), 'LineWidth', 2);
        legendLabel = cell(params.n_mts + 1, 1); %, 'Crosslinkers', 'Combined'};
        for i_pf = 1 : 1 : params.n_mts
           legendLabel{i_pf} = sprintf('Protofilament %i', int32(i_pf)); 
        end
        legendLabel{params.n_mts + 1} = 'Average across all';
        legend(legendLabel, 'Location', 'northeastoutside');

        dim = [0.7425 0.3 .3 .3];
        time = i * params.time_per_datapoint;
        str = sprintf('Time: %i seconds', int32(time));
        annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on');

        drawnow();
        frame = getframe(gcf);
        writeVideo(v, getframe(gcf));
    end

end

close(v);
