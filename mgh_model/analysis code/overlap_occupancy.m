clear all

% Parameters from sim
n_datapoints = 100000;
motor_ID = 2;
mt_length = 500;
n_mts = 2;
xlink_cutoff = 7;

% File info
simName = 'presS4';
fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
mtFileName = '%s_MTcoord.file';
motorFileName = '%s_motorID.file';
xlinkFileName = '%s_xlinkID.file';
tethFileName = '%s_tether_coord.file';
mtFile = sprintf(fileDirectory, sprintf(mtFileName, simName));
motorFile = sprintf(fileDirectory, sprintf(motorFileName, simName));
xlinkFile = sprintf(fileDirectory, sprintf(xlinkFileName, simName));
tethFile = sprintf(fileDirectory, sprintf(tethFileName, simName));

% Figure parameters (i.e., how they appear)
n_frames = 100000;
frames_per_plot = 100;
start_frame = 001;
site_height = 1;
site_width = 1;

% Videowriter details
v = VideoWriter('overlap_expansion.avi');
v.FrameRate = (n_frames / frames_per_plot) / 25;
open(v);
frame_box = [0 0 1545 200];

% Figure details
fig1 = figure;
set(fig1, 'Position', [0 100 1600 400])

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

% Run through all datapoints; each one is a frame in our movie
for i_data=start_frame:frames_per_plot:(start_frame + n_frames - 1)
    
    bot_MT_coord = mt_data(:, 1, i_data);
    top_MT_coord = mt_data(:, 2, i_data);
    overlap_start_coord;
    overlap_end_coord;
    if(bot_MT_coord < top_mt_coord)
        overlap_start_coord = top_mt_coord;
        overlap_end_coord = bot_mt_coord + mt_length;
    else
        overlap_start_coord = bot_mt_coord;
        overlap_end_coord = top_mt_coord + mt_length;
    end
    
    overlap_length = overlap_end_coord - overlap_end_coord;
    
    n_overlap_occupancy = 0; 
    
    for i_mt=1:1:n_mts
        motor_IDs = motor_data(:, i_mt, i_data);
        MT_coord = mt_data(:, i_mt, i_datapoint);
        for i_motor=1:1:mt_length - 1
            motor_coord = MT_coord + i_motor;
            if(motor_coord <= overlap_start_coord && motor_coord >= overlap_end_coord)
                if(motor_IDs(i_motor) ~= -1 && motor_IDs(i_motor) == motor_IDs(i_motor + 1))
                    n_overlap_occupancy = n_overlap_occupancy + 1;   
                end
            end
        end   
    end
    %{
        xlink_IDs = xlink_data(:, i_mt, i_data);
        neighb_xlink_IDs;
        if(i_mt == 1)
            ne
            
        end
        
      %}
end