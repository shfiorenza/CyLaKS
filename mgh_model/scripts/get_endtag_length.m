function endtag_length = get_endtag_length(sim_name, n_sites)

motor_speciesID = 2;
xlink_speciesID = 1;
n_datapoints = 10000;
starting_point = 5000;
active_datapoints = n_datapoints - starting_point;
fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
fileStruct = '%s_%i_occupancy.file';

fileName = sprintf(fileDirectory, sprintf(fileStruct, sim_name, n_sites));
data_file = fopen(fileName);
motor_raw_data = fread(data_file, [n_sites, n_datapoints], '*int');
xlink_raw_data = motor_raw_data;
fclose(data_file);

motor_raw_data(motor_raw_data ~= motor_speciesID) = 0;
motor_raw_data(motor_raw_data == motor_speciesID) = 1;
xlink_raw_data(xlink_raw_data ~= xlink_speciesID) = 0;
xlink_raw_data(xlink_raw_data == xlink_speciesID) = 1;

motor_avg_occupancy = zeros([n_sites 1]);
xlink_avg_occupancy = zeros([n_sites 1]);

% Read in and average occupancy data over all datapoints
for i=starting_point:1:n_datapoints
    motor_avg_occupancy(:,1) = motor_avg_occupancy(:,1) + double(motor_raw_data(:,i))./active_datapoints;
    xlink_avg_occupancy(:,1) = xlink_avg_occupancy(:,1) + double(xlink_raw_data(:,i))./active_datapoints;
end

smooth_window = n_sites / 10;
motor_occupancy = smoothdata(motor_avg_occupancy, 'movmean', smooth_window);
xlink_occupancy = smoothdata(xlink_avg_occupancy, 'movmean', smooth_window);
net_occupancy = motor_occupancy + xlink_occupancy;
occupancy_slope = smoothdata(gradient(net_occupancy, 0.008),'movmean', smooth_window);
occupancy_accel = smoothdata(gradient(occupancy_slope, 0.008), 'movmean', smooth_window);
max_occupancy = max(net_occupancy);
min_slope = min(occupancy_slope);

past_threshold = false;
i_threshold = 0;
endtag_site = 0;
for i_site=1:n_sites
    if(~past_threshold && net_occupancy(i_site) < 0.5*max_occupancy)
        past_threshold = true;
        i_threshold = i_site;
    end
    if(past_threshold && occupancy_slope(i_site) > (min_slope / 2))
        endtag_site = i_site;
        break;
    end
end
if endtag_site == 0
    [~, i_peak] = max(occupancy_accel(i_threshold:n_sites));
    endtag_site = i_threshold + i_peak;
end
endtag_length = endtag_site*0.008;

end
