
clear all;
% Often-changed variables
n_sites = 50000;
runlengths_used = [1.2, 6.0, 12.0, 18.0];
suffixes = ["_50","","_170"];
n_suffixes = length(suffixes);
% Pseudo-constant variables
n_mts = 1;
delta_t = 0.0001;
n_steps = 10000000;
n_datapoints = 10000;
time_per_datapoint = delta_t * n_steps / n_datapoints;
starting_point = 1;
active_datapoints = n_datapoints - starting_point;
time_cutoff = 0; %0.3;    %in seconds
site_size = 0.008; % in um

n_runlengths = length(runlengths_used);
expected_runlengths = zeros(n_runlengths, n_suffixes);
for i_suffix=1:1:n_suffixes
    expected_runlengths(:, i_suffix) = runlengths_used;
end
actual_runlengths = zeros(n_runlengths, n_suffixes);
runlength_errors = zeros(n_runlengths, n_suffixes);
actual_velocities = zeros(n_runlengths, n_suffixes);
velocity_errors = zeros(n_runlengths, n_suffixes);

fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/output_pro/%s';
motorFileStruct = '%s_motorID.file';

for i_suffix=1:1:n_suffixes
    
    for i_runlength=1:1:n_runlengths
        runlength = expected_runlengths(i_runlength, i_suffix);
        simName = sprintf('Processivity_%#.1f%s', runlength, suffixes(i_suffix));
        
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
                            active_motors(i_mt, n_active(i_mt)) = motor_ID;
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
        
        % matlab's exponential fit always goes to zero; offset it appropriately
        min_run = min(run_lengths);
        run_lengths = run_lengths - min_run;
        % plot run length histogram
        n_bins = int32(sqrt(n_runs));
        %   hist = histfit(run_lengths, n_bins, 'exponential');
        % fit distribution of shifted run lengths
        run_dist = fitdist(run_lengths, 'exponential');
        mean_run = run_dist.mu + min_run;
        conf_inv_run = paramci(run_dist);
        sigma_run = abs(conf_inv_run(2) - conf_inv_run(1)) / 2;
        
        actual_runlengths(i_runlength, i_suffix) =  mean_run;
        runlength_errors(i_runlength, i_suffix) = sigma_run;
        
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
        
        vel = round(single(mean_run * 1000 / mean_time), -1);
        sigma_vel = round(single(sqrt((sigma_time/mean_time)^2 + (sigma_run/mean_run)^2) * vel),-1);
        
        actual_velocities(i_runlength, i_suffix) = vel;
        velocity_errors(i_runlength, i_suffix) = sigma_vel;
    end
end
%plot(theoretical_runlengths, actual_runlengths, 'LineWidth', 2);
fig1 = figure(1);
set(fig1,'Position', [50, 50, 1.5*600, 1.5*500])
plot_one = subplot(2, 1, 1);
set(plot_one, 'Position', [.1 .55 .8 .4]);
errorbar(expected_runlengths, actual_runlengths, runlength_errors, ...
    'o','LineWidth',2);
hold on
plot(expected_runlengths, expected_runlengths, '--r', 'LineWidth', 1);
% Cosmetic stuff
title('Processivity scaling vs linear theory');
%xlabel('Theoretical run length (um)');
ylabel('Measured run length (um)');
xlim([0 expected_runlengths(n_runlengths)+2]);
ylim([0 expected_runlengths(n_runlengths)+2]);
xticks([0,5,10,15,20]);
yticks([0,5,10,15,20]);
legend({'k\_hydrolyze = 50 s^{-1}','k\_hydrolyze = 110 s^{-1}', ... 
    'k\_hydrolyze = 170 s^{-1}', 'Predicted from theory'}, 'location', 'southeast');

plot_two = subplot(2, 1, 2);
set(plot_two, 'Position', [.1 .08 .8 .4]);
errorbar(expected_runlengths, actual_velocities, velocity_errors, ...
    'o-','LineWidth',2);
%title('Processivity scaling vs linear theory');
xlabel('Set run length (um)');
ylabel('Measured stepping velocity (nm/s)');
xlim([0 expected_runlengths(n_runlengths)+2]);
ylim([0 1200]);
xticks([0,5,10,15,20]);
yticks([0,300,600,900,1200]);