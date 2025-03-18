function endtag_length = get_endtag_length(sim_name)
    plot_flag = true;

    motor_speciesID = 2;
    xlink_speciesID = 1;
    site_size = 0.0082; % in um
    % Load parameter structure
    params = load_parameters(sim_name);

    occupancy_filename = sprintf('%s_occupancy.file', sim_name);
    occupancy = zeros(params.max_sites, params.n_mts, params.n_datapoints) - 1;
    occupancy = load_data(occupancy, occupancy_filename, '*int');
    xlink_raw_data = occupancy;
    motor_raw_data = occupancy;

    motor_avg_occupancy = zeros([params.max_sites 1]);
    motor_raw_data(motor_raw_data ~= motor_speciesID) = 0;
    motor_raw_data(motor_raw_data == motor_speciesID) = 1;

    xlink_avg_occupancy = zeros([params.max_sites 1]);
    xlink_raw_data(xlink_raw_data ~= xlink_speciesID) = 0;
    xlink_raw_data(xlink_raw_data == xlink_speciesID) = 1;

    dwell_time = 100;
    dwell_steps = int32(dwell_time / params.time_per_datapoint);

    starting_point = params.n_datapoints - dwell_steps; % 1
    active_datapoints = double(dwell_steps); % n_datapoints - starting_point + 1;

    % Read in and average occupancy data over all datapoints
    for i = starting_point:1:params.n_datapoints
        for i_mt = 1 : params.n_mts
            motor_avg_occupancy(:, 1) = motor_avg_occupancy(:, 1) + double(motor_raw_data(:, i_mt, i)) ./ (active_datapoints * params.n_mts);
            xlink_avg_occupancy(:, 1) = xlink_avg_occupancy(:, 1) + double(xlink_raw_data(:, i_mt, i)) ./ (active_datapoints * params.n_mts);
        end
    end

    smooth_window = 32; % should be equivalent to diffraction limit%n_sites / 20;
    motor_occupancy = smoothdata(motor_avg_occupancy, 'movmean', smooth_window);
    xlink_occupancy = smoothdata(xlink_avg_occupancy, 'movmean', smooth_window);
    %net_occupancy = motor_occupancy + xlink_occupancy;
    net_occupancy = xlink_occupancy;
    occupancy_slope = smoothdata(gradient(net_occupancy, 0.008), 'movmean', smooth_window);
    occupancy_accel = smoothdata(gradient(occupancy_slope, 0.008), 'movmean', smooth_window);
    max_occupancy = max(net_occupancy);
    min_slope = min(occupancy_slope);

    past_threshold = false;
    i_threshold = 0;
    endtag_site = 0;

    for i_site = 1:params.max_sites 
        if (~past_threshold && net_occupancy(i_site) < 0.9 * max_occupancy)
            past_threshold = true;
            i_threshold = i_site;
        end

        if (past_threshold && occupancy_slope(i_site) > (0.5 * min_slope))
            endtag_site = i_site;
            break;
        end
        %}
    end
    %{
    if endtag_site == 0
        [max_accel, i_peak] = max(occupancy_accel(i_threshold + 1:params.max_sites));
        endtag_site = i_threshold + i_peak;
        disp('boop')
    end
    %}
    if plot_flag == true
        fig = figure();
        plot(net_occupancy);
        hold on
        plot([endtag_site endtag_site], [-1 2], '--')
        ylim([0 1])
        title(sim_name)
        drawnow
    end
    endtag_length = endtag_site * site_size;

end
