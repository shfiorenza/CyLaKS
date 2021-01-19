function endtag_length = get_endtag_length(sim_name)

    %sim_name = "/home/shane/Projects/overlap_analysis/mgh_model/endtag_250_0";

    motor_speciesID = 2;
    xlink_speciesID = 1;
    site_size = 0.0082; % in um

    % Open log file and parse it into param labels & their values
    log_file = sprintf('%s.log', sim_name);
    log = textscan(fileread(log_file), '%s %s', 'Delimiter', '=');
    params = log{1, 1};
    values = log{1, 2};
    n_sites = values{contains(params, "n_sites[0]")};
    n_sites = sscanf(n_sites, '%i');
    % Read in system params
    dt = sscanf(values{contains(params, "dt ")}, '%g');
    steps_per_datapoint = str2double(values{contains(params, "n_steps_per_snapshot ")});
    time_per_datapoint = dt * steps_per_datapoint;
    n_datapoints = str2double(values{contains(params, "n_datapoints ")});
    % Use actual recorded number of datapoints to parse thru data/etc
    if any(contains(params, "N_DATAPOINTS ") ~= 0)
        n_datapoints = str2double(values{contains(params, "N_DATAPOINTS ")});
    end
   
    fileName = sprintf("%s_occupancy.file", sim_name);
    data_file = fopen(fileName);
    motor_raw_data = fread(data_file, [n_sites, n_datapoints], '*int');
    xlink_raw_data = motor_raw_data;
    fclose(data_file);

    motor_avg_occupancy = zeros([n_sites 1]);
    motor_raw_data(motor_raw_data ~= motor_speciesID) = 0;
    motor_raw_data(motor_raw_data == motor_speciesID) = 1;

    xlink_avg_occupancy = zeros([n_sites 1]);
    xlink_raw_data(xlink_raw_data ~= xlink_speciesID) = 0;
    xlink_raw_data(xlink_raw_data == xlink_speciesID) = 1;

    dwell_time = 10;
    dwell_steps = int32(dwell_time / time_per_datapoint);

    starting_point = 1; % n_datapoints - dwell_steps;% 1
    active_datapoints = n_datapoints; %double(dwell_steps);% n_datapoints - starting_point + 1;

    % Read in and average occupancy data over all datapoints
    for i = starting_point:1:n_datapoints
        motor_avg_occupancy(:, 1) = motor_avg_occupancy(:, 1) + double(motor_raw_data(:, i)) ./ active_datapoints;
        xlink_avg_occupancy(:, 1) = xlink_avg_occupancy(:, 1) + double(xlink_raw_data(:, i)) ./ active_datapoints;
    end

    smooth_window = 32; % should be equivalent to diffraction limit%n_sites / 20;
    motor_occupancy = smoothdata(motor_avg_occupancy, 'movmean', smooth_window);
    xlink_occupancy = smoothdata(xlink_avg_occupancy, 'movmean', smooth_window);
    net_occupancy = motor_occupancy + xlink_occupancy;
    occupancy_slope = smoothdata(gradient(net_occupancy, 0.008), 'movmean', smooth_window);
    occupancy_accel = smoothdata(gradient(occupancy_slope, 0.008), 'movmean', smooth_window);
    max_occupancy = max(net_occupancy);
    min_slope = min(occupancy_slope);

    past_threshold = false;
    i_threshold = 0;
    endtag_site = 0;

    for i_site = 1:n_sites

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
        [max_accel, i_peak] = max(occupancy_accel(i_threshold + 1:n_sites));
        endtag_site = i_threshold + i_peak;
    end

    endtag_length = endtag_site * site_size;

end
