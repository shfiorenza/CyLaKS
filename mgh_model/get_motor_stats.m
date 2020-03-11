function mot_stats = get_motor_stats(sim_name)

% Often-changed variables
n_sites = 5000;
simName = sim_name;
% Pseudo-constant variables
n_mts = 1;
delta_t = 0.000025;
n_steps = 40000000;
n_datapoints = 100000;
time_per_datapoint = delta_t * n_steps / n_datapoints;
starting_point = 1;
time_cutoff = time_per_datapoint; %0.3;    %in seconds
site_size = 0.008; % in um

fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
motorFileStruct = '%s_motorID.file';

motorFileName = sprintf(fileDirectory, sprintf(motorFileStruct, simName));
motor_data_file = fopen(motorFileName);
raw_motor_data = fread(motor_data_file, n_mts * n_sites * n_datapoints, '*int');
fclose(motor_data_file);
motor_data = reshape(raw_motor_data, n_sites, n_mts, n_datapoints);

% have an active list for each MT
active_motors = zeros([n_mts n_mts*n_sites]);
n_active = zeros([n_mts 1]);

% motor ID is unique; make following arrays serial w/ ID as index
runlengths = zeros([(n_mts*n_sites) 1]);
lifetimes = zeros([(n_mts*n_sites) 1]);
velocities = zeros([(n_mts*n_sites) 1]);
n_runs = 0;
starting_site = zeros([n_mts*n_sites 1]) - 1;
starting_datapoint = zeros([n_mts*n_sites 1]) - 1;

for i_data = starting_point:1:n_datapoints - 1
    for i_mt = 1:1:n_mts
        motor_IDs = motor_data(:, i_mt, i_data);
        future_IDs = motor_data(:, i_mt, i_data + 1);
        endtag_boundary = 1;
        % Determine end-tag region; ignore motors that terminate here
        for i_site=1:n_sites
           motor_ID = future_IDs(i_site);
           if motor_ID ~= -1
               endtag_boundary = i_site + 1;
           else
               break;    
           end
        end
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
                % Otherwise if a motor is found, only count first head
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

        n_deleted = 0;
        for i_motor = 1:1:n_active(i_mt)
            i_adj = i_motor - n_deleted;
            motor_ID = active_motors(i_mt, i_adj);
            future_site = find(future_IDs == motor_ID, 1);
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
                velocity = (run_length / run_time) * 1000; % convert to nm/s
                % If time bound is above time cutoff, add to data
                if run_time > time_cutoff && end_site(1) > endtag_boundary
                    n_runs = n_runs + 1;
                    runlengths(n_runs) = run_length;
                    lifetimes(n_runs) = run_time;
                    velocities(n_runs) = velocity;  
                end
                starting_site(motor_ID) = -1;
                starting_datapoint(motor_ID) = -1;
                % Switch now-deleted entry with last entry in active_motors
                active_motors(i_mt, i_adj) = active_motors(i_mt, n_active(i_mt));
                active_motors(i_mt, n_active(i_mt)) = -1;
                n_active(i_mt) = n_active(i_mt) - 1;
                n_deleted = n_deleted + 1;
            end
        end
    end
end
% trim arrays to get rid of un-used containers
runlengths = runlengths(1:n_runs);
lifetimes = lifetimes(1:n_runs);
velocities = velocities(1:n_runs);

% matlab's exponential fit always goes to zero; offset it appropriately
min_run = min(runlengths);
runlengths = runlengths - min_run;
% fit distribution of shifted run lengths
run_dist = fitdist(runlengths, 'exponential');
mean_run = run_dist.mu + min_run;

% matlab's exponential fit always goes to zero; offset it appropriately
min_time = min(lifetimes);
lifetimes = lifetimes - min_time;
% fit distribution of shifted run times
time_dist = fitdist(lifetimes, 'exponential');
% get mean run time
mean_time = time_dist.mu + min_time;

vel_dist = fitdist(velocities, 'normal');
mean_vel = vel_dist.mu;


mot_stats(1) = mean_run;
mot_stats(2) = mean_time;
mot_stats(3) = mean_vel;

end
