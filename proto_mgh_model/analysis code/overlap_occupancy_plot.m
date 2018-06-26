clear all

% Parameters from sim
n_datapoints = 100000;
%motor_ID = 2;
mt_length = 1000;
n_mts = 2;
xlink_cutoff = 7;

% File info
simName = 'testlong';
fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
mtFileName = '%s_mt_coord.file';
motorFileName = '%s_motorID.file';
xlinkFileName = '%s_xlinkID.file';
tethFileName = '%s_tether_coord.file';
mtFile = sprintf(fileDirectory, sprintf(mtFileName, simName));
motorFile = sprintf(fileDirectory, sprintf(motorFileName, simName));
xlinkFile = sprintf(fileDirectory, sprintf(xlinkFileName, simName));
tethFile = sprintf(fileDirectory, sprintf(tethFileName, simName));

mt_data_file = fopen(mtFile);
mt_raw_data = fread(mt_data_file, [n_mts * n_datapoints], '*double');
fclose(mt_data_file);
mt_data = reshape(mt_raw_data, n_mts, n_datapoints);

motor_data_file = fopen(motorFile);
motor_raw_data = fread(motor_data_file, [n_mts * mt_length * n_datapoints], '*int');
fclose(motor_data_file);
motor_data = reshape(motor_raw_data, mt_length, n_mts, n_datapoints);

xlink_data_file = fopen(xlinkFile);
xlink_raw_data = fread(xlink_data_file, [n_mts * mt_length * n_datapoints], '*int');
fclose(xlink_data_file);
xlink_data = reshape(xlink_raw_data, mt_length, n_mts, n_datapoints);

overlap_occupancy_data = zeros(n_datapoints, 1);
overlap_length_data = zeros(n_datapoints, 1);

% Run through all datapoints, record overlap length 
% and occupancy of force-exerting proteins
for i_data=1:1:n_datapoints
    
    bot_mt_coord = mt_data(1, i_data);
    top_mt_coord = mt_data(2, i_data);
    %overlap_start_coord = 0;
    %overlap_end_coord = 0;
    if(bot_mt_coord < top_mt_coord)
        overlap_start_coord = top_mt_coord;
        overlap_end_coord = bot_mt_coord + mt_length;
    else
        overlap_start_coord = bot_mt_coord;
        overlap_end_coord = top_mt_coord + mt_length;
    end
    
    overlap_length = overlap_end_coord - overlap_start_coord;
    
    overlap_occupancy = 0; 
    
    for i_mt=1:1:n_mts
        motor_IDs = motor_data(:, i_mt, i_data);
        xlink_IDs = xlink_data(:, i_mt, i_data);
        MT_coord = mt_data(i_mt, i_data);
        for i_motor=1:1:mt_length - 1
            motor_coord = MT_coord + i_motor;
            if(motor_coord < overlap_end_coord & motor_coord > overlap_start_coord)
                if(motor_IDs(i_motor) ~= -1 & motor_IDs(i_motor) == motor_IDs(i_motor + 1))
                    overlap_occupancy = overlap_occupancy + 1;   
                end
            end
        end   
        for i_xlink=1:1:mt_length
            xlink_coord = MT_coord + i_xlink;
            if(xlink_coord < overlap_end_coord & xlink_coord > overlap_start_coord)
                if(xlink_IDs(i_xlink ~= -1))
                    overlap_occupancy = overlap_occupancy + 1;
                end
            end
        end
    end
    
    overlap_occupancy_data(i_data) = overlap_occupancy;
    overlap_length_data(i_data) = overlap_length;
    
end

figure(1);
plot(linspace(0, 1, n_datapoints), overlap_occupancy_data, 'Color', 'b');
hold on
plot(linspace(0, 1, n_datapoints), overlap_length_data, 'Color', 'r');
legend({'Overlap occupancy', 'Overlap length'}, 'Location', 'northeastoutside');
ylim([0 inf]);