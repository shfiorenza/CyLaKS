
clear all;
% Often-changed variables
n_sites = 5000;
simName = 'MotorAppForce_0.0pN';
% Pseudo-constant variables
n_mts = 1;
delta_t = 0.0001;
n_steps = 1000000; %0
n_datapoints = 10000;
time_per_datapoint = delta_t * n_steps / n_datapoints;
starting_point = 1;
active_datapoints = n_datapoints - starting_point;
time_cutoff = 0; %0.3;    %in seconds
site_size = 0.008; % in um

fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
motorFileStruct = '%s_motorID.file';

f_applied = [0, 1, 2, 3, 4, 5, 6];
mean_run_lengths = zeros([length(f_applied) 1]);
mean_velocities = zeros([length(f_applied) 1]);

for i_force=1:1:length(f_applied)
    
    simName = sprintf('MotorAppForce_%#.2gpN_b', f_applied(i_force));
    
    motorFileName = sprintf(fileDirectory, sprintf(motorFileStruct, simName)); 
    motor_data_file = fopen(motorFileName);
    raw_motor_data = fread(motor_data_file, [n_mts * n_sites * n_datapoints], '*int');
    fclose(motor_data_file);
    motor_data = reshape(raw_motor_data, n_sites, n_mts, n_datapoints);
    
    % have an active list for each MT
    active_motors = zeros([n_mts n_mts*n_sites]);
    n_active = zeros([n_mts 1]);
    
    % motor ID is unique; make following arrays serial w/ ID as index
    run_lengths = zeros([(n_mts*n_sites*n_datapoints) 1]);
    run_times = zeros([(n_mts*n_sites*n_datapoints) 1]);
    n_runs = 0;
    starting_site = zeros([n_mts*n_sites 1]) - 1;
    starting_datapoint = zeros([n_mts*n_sites 1]) - 1;
    
    for i_data = starting_point:1:n_datapoints - 1
        for i_mt = 1:1:n_mts
            motor_IDs = motor_data(:, i_mt, i_data);
            % Scan through IDs of bound motors (-1 means no motor on that site)
            for i_site = 1:1:n_sites
                motor_ID = motor_IDs(i_site);
                % Always count motor on first site
                if motor_ID > 0 && i_site == 1
                    % Record the motor's starting site if this is the first time
                    % seeing it (-1 means it was not seen last datapoint)
                    if starting_site(motor_ID) == -1
                        starting_site(motor_ID) = i_site;
                        starting_datapoint(motor_ID) = i_data;
                        n_active(i_mt) = n_active(i_mt) + 1;
                        active_motors(i_mt, n_active) = motor_ID;
                    end
                    % Otherwise if a motor is found, only count one head
                elseif motor_ID > 0 && motor_IDs(i_site - 1) ~= motor_ID
                    % Record the motor's starting site if this is the first time
                    % seeing it (-1 means it was not seen last datapoint)
                    if starting_site(motor_ID) == -1
                        starting_site(motor_ID) = i_site;
                        starting_datapoint(motor_ID) = i_data;
                        n_active(i_mt) = n_active(i_mt) + 1;
                        active_motors(i_mt, n_active(i_mt)) = motor_ID;
                    end
                end
            end
            % Check one datapoint into the future to see if any motors unbound
            future_IDs = motor_data(:, i_mt, i_data + 1);
            n_deleted = 0;
            for i_motor = 1:1:n_active(i_mt)
                i_adj = i_motor - n_deleted;
                motor_ID = active_motors(i_mt, i_adj);
                future_site = find(future_IDs == motor_ID);
                if isempty(future_site)
                    % Calculate distance traveled
                    end_site = find(motor_IDs == motor_ID);
                    start_site = starting_site(motor_ID);
                    delta = abs(end_site(1) - start_site);
                    run_length = delta * site_size;
                    % Calculate time bound
                    start_datapoint = starting_datapoint(motor_ID);
                    delta_t = abs(i_data - start_datapoint);
                    run_time = delta_t * time_per_datapoint;
                    % If time bound is above time cutoff, add to data
                    if run_time >= time_cutoff
                        n_runs = n_runs + 1;
                        run_lengths(n_runs) = run_length;
                        run_times(n_runs) = run_time;
                    end            
                    starting_site(motor_ID) = -1;
                    starting_datapoint(motor_ID) = -1;
                    active_motors(i_mt, i_adj) = active_motors(i_mt, n_active(i_mt));
                    active_motors(i_mt, n_active(i_mt)) = -1;
                    n_active(i_mt) = n_active(i_mt) - 1;
                    n_deleted = n_deleted + 1;
                end
            end
        end
    end
    % trim arrays to get rid of un-used containers
    run_lengths = run_lengths(1:n_runs);
    run_times = run_times(1:n_runs);
    % prep figure   
%    fig1 = figure();
%    set(fig1, 'Position', [50, 50, 960, 600]);
    
    % matlab's exponential fit always goes to zero; offset it appropriately
    min_run = min(run_lengths);
    run_lengths = run_lengths - min_run;
    % plot run length histogram
    n_bins = int32(sqrt(n_runs));
%    hist = histfit(run_lengths, n_bins, 'exponential');  
    % fit distribution of shifted run lengths
    run_dist = fitdist(run_lengths, 'exponential');
    mean_run = run_dist.mu + min_run;
    conf_inv_run = paramci(run_dist);
    sigma_run = abs(conf_inv_run(2) - conf_inv_run(1)) / 2;
    
    % matlab's exponential fit always goes to zero; offset it appropriately
    min_time = min(run_times);
    run_times = run_times - min_run;
    % fit distribution of shifted run times
    time_dist = fitdist(run_times, 'exponential');
    % get mean run time
    mean_time = time_dist.mu;
    % get confidence interval of mean run time
    conf_inv_time = paramci(time_dist); 
    % calculate sigma
    sigma_time = abs(conf_inv_time(2) - conf_inv_time(1)) / 2;
    
    vel = mean_run * 1000 / mean_time;
    %{
    % Display mean values on run length histogram itself
    dim1 = [0.55 0.65 0.2 0.2];
    str1 = sprintf('Mean run length: %#.3g +/- %#.1g microns', mean_run, sigma_run);
    annotation('textbox',dim1,'String',str1,'FitBoxToText','on');
    dim2 = [0.55 0.6 0.2 0.2];
    str2 = sprintf('Mean run time: %#.3g +/- %#.1g seconds', mean_time, sigma_time);
    annotation('textbox',dim2,'String',str2,'FitBoxToText','on'); 
    dim3 = [0.55 0.55 0.2 0.2];
    vel = round(single(mean_run * 1000 / mean_time), -1);
    sigma_vel = round(single(sqrt((sigma_time/mean_time)^2 + (sigma_run/mean_run)^2) * vel),-1);
    str3 = sprintf('Calculated velocity: %#d +/- %#d nm/s', vel, sigma_vel);
    annotation('textbox',dim3,'String',str3,'FitBoxToText','on'); 
    % Cosmetic stuff
    title(sprintf('Run length histogram for %i micron MT', int32(n_sites * 0.008))); 
    xlabel('Run length (um)');
    ylabel('Counts');
    %}
    
    mean_run_lengths(i_force) = mean_run;
    mean_velocities(i_force) = vel;
end

scaled_runs_expected = [1.2,0.76,0.44,0.273,0.164,0.098,0.0545];
scaled_vel_expected = [640,597,555,469,341,171,100];

fig1 = figure();
set(fig1, 'Position', [50, 50, 960, 600])
sgtitle('unbinding same as paper (sigma x2.5); v ~ (1 - F_x/F_s)^{1/5}');
subplot(1, 2, 1);
plot(f_applied, mean_run_lengths, 'LineWidth', 2);
hold on
plot(f_applied, scaled_runs_expected, 'LineWidth', 2);
title('Processivity scaling with applied force'); 
xlabel('Force applied (pN)');
ylabel('Mean processivity (um)');
grid on
grid minor
legend({'Simulation', 'Expected'}, 'location', 'northeast');
subplot(1, 2, 2);
plot(f_applied, mean_velocities, 'LineWidth', 2);
hold on
plot(f_applied, scaled_vel_expected, 'LineWidth', 2);
title('Velocity scaling with applied force'); 
xlabel('Force applied (pN)');
ylabel('Mean stepping velocity (nm/s)');
grid on
grid minor
legend({'Simulation', 'Expected'}, 'location', 'northeast');