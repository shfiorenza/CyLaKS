clear variables;
% Often-changed variables

sim_name = 'hybrid_motor_100_0';
steps_per_plot = 10;
movie_name = 'test2';
file_dir = '/home/shane/projects/CyLaKS';
%file_dir='.';

% Open log file and parse it into param labels & their values
log_file = sprintf('%s/%s', file_dir, sprintf('%s.log', sim_name));
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


% Pseudo-constant variables
motor_speciesID = 2;
xlink_speciesID = 1;

fileStruct = '%s_occupancy.file';
legendLabel = {'Fractional occupancy of motors'}; %, 'Crosslinkers', 'Combined'};
% Videowriter details
v = VideoWriter(movie_name);
v.FrameRate = (n_datapoints / steps_per_plot) / 15;
open(v);
frame_box = [0, 0, 1.5 * 480, 1.5 * 300];

fileName = sprintf("%s/%s", file_dir, sprintf(fileStruct, sim_name));
data_file = fopen(fileName);
motor_raw_data = fread(data_file, [max_sites, n_datapoints], '*int');
xlink_raw_data = motor_raw_data;
fclose(data_file);

motor_avg_occupancy = zeros([max_sites 1]);
xlink_avg_occupancy = zeros([max_sites 1]);

motor_raw_data(motor_raw_data ~= motor_speciesID) = 0;
motor_raw_data(motor_raw_data == motor_speciesID) = 1;
xlink_raw_data(xlink_raw_data ~= xlink_speciesID) = 0;
xlink_raw_data(xlink_raw_data == xlink_speciesID) = 1;

fig1 = figure(1);
set(fig1, 'Position', [50, 50, 1200, 600])

% Read in and average occupancy data over all datapoints
for i = 1:1:int32(n_datapoints)
    motor_avg_occupancy(:, 1) = motor_avg_occupancy(:, 1) + double(motor_raw_data(:, i)) ./ steps_per_plot;
    xlink_avg_occupancy(:, 1) = xlink_avg_occupancy(:, 1) + double(xlink_raw_data(:, i)) ./ steps_per_plot;

    if (mod(i, steps_per_plot) == 0)
        smooth_window = 32; % should be equivalent to diffraction limit
        motor_occupancy = smoothdata(motor_avg_occupancy, 'movmean', smooth_window);
        xlink_occupancy = smoothdata(xlink_avg_occupancy, 'movmean', smooth_window);
        net_occupancy = motor_occupancy + xlink_occupancy;
        occupancy_slope = smoothdata(gradient(net_occupancy, 0.008), 'movmean', smooth_window);
        occupancy_accel = smoothdata(gradient(occupancy_slope, 0.008), 'movmean', smooth_window);
        max_occupancy = max(net_occupancy);
        min_slope = min(occupancy_slope);

        motor_avg_occupancy(:, 1) = 0;
        xlink_avg_occupancy(:, 1) = 0;

        past_threshold = false;
        i_threshold = 0;
        endtag_site = 0;

        for i_site = 1:max_sites

            if (~past_threshold && net_occupancy(i_site) < 0.5 * max_occupancy)
                past_threshold = true;
                i_threshold = i_site;
            end

            if (past_threshold && occupancy_slope(i_site) > (0.5 * min_slope))
                endtag_site = i_site;
                break;
            end

        end

        if endtag_site == 0
            [max_accel, i_peak] = max(occupancy_accel(i_threshold + 1:max_sites));
            endtag_site = i_threshold + i_peak;
        end

        endtag_length = endtag_site * site_size;

        %%plot fig%%
        clf;
        ax = axes('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
        hold all

        plot(linspace(0, max_sites * site_size, max_sites), motor_occupancy);
        %plot(linspace(0, n_sites*0.008, n_sites), xlink_occupancy);
        %plot(linspace(0, n_sites*0.008, n_sites), net_occupancy);
        % plot(linspace(0, n_sites*0.008, n_sites), occupancy_slope);
        % plot(linspace(0, n_sites*0.008, n_sites), occupancy_accel);
        % Put vertical red line where endtag starting position is
        plot([endtag_length endtag_length], [0 1], ':', 'LineWidth', 2, 'Color', 'red');

        title(sprintf('End-tag length: %g microns', endtag_length), 'FontSize', 14)
        xlabel('Distance from plus-end (\mum)', 'FontSize', 14);
        ylabel('Fractional site occupancy', 'FontSize', 14);
        ylim([0 1]);
        set(gca, 'FontSize', 14);
        axis = gca;
        %axis.TickDir = 'out';
        axis.Box = 'off';
        axis.GridLineStyle = '-';
        set(findall(axis, 'Type', 'Line'), 'LineWidth', 2);
        %legend(legendLabel, 'Location', 'northeast');

        dim = [0.75 0.55 .3 .3];
        time = i * time_per_datapoint;
        str = sprintf('Time: %i seconds', int32(time));
        annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on');

        drawnow();
        frame = getframe(gcf);
        writeVideo(v, getframe(gcf));
    end

end

close(v);
