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

motor_avg_occupancy = zeros([n_sites 1]);
xlink_avg_occupancy = zeros([n_sites 1]);

motor_raw_data(motor_raw_data ~= motor_speciesID) = 0;
motor_raw_data(motor_raw_data == motor_speciesID) = 1;

xlink_raw_data(xlink_raw_data ~= xlink_speciesID) = 0;
xlink_raw_data(xlink_raw_data == xlink_speciesID) = 1;

% Read in and average occupancy data over all datapoints
for i=starting_point:1:n_datapoints
    motor_avg_occupancy(:,1) = motor_avg_occupancy(:,1) + double(motor_raw_data(:,i))./active_datapoints;
    xlink_avg_occupancy(:,1) = xlink_avg_occupancy(:,1) + double(xlink_raw_data(:,i))./active_datapoints;
end

motor_smoothed_avg = smoothdata(motor_avg_occupancy);
xlink_smoothed_avg = smoothdata(xlink_avg_occupancy);
net_smoothed_avg = motor_smoothed_avg + xlink_smoothed_avg;
occupancy_slope = abs(gradient(net_smoothed_avg, 0.08));
max_slope = max(occupancy_slope);

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

endtag_length = double(endtag_site*0.008);

end
