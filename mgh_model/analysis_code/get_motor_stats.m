function mot_stats = get_motor_stats(sim_name)

% Open log file and parse it into param labels & their values
log_file = sprintf('%s.log', sim_name);
log = textscan(fileread(log_file),'%s %s', 'Delimiter', '=');
params = log{1,1};
values = log{1,2};
% Read in number of MTs
n_mts = str2double(values{contains(params, 'count')});
n_sites = values{contains(params, 'length')};
n_sites = sscanf(n_sites, '%i');
% Read in system params
delta_t = sscanf(values{contains(params, 'delta_t')}, '%g');
total_steps = str2double(values{contains(params, 'n_steps')});
data_threshold = sscanf(values{contains(params, 'data_threshold')}, '%g');
if any(contains(params, 'DATA_THRESHOLD') ~= 0)
   data_threshold = str2double(values{contains(params, 'DATA_THRESHOLD')});
end
n_steps = total_steps - data_threshold;
% Use max possible number of datapoints to calculate time_per_datapoint (as is done in Sim)
n_datapoints = str2double(values{contains(params, 'n_datapoints')});
time_per_datapoint = delta_t * n_steps / n_datapoints;
site_size = 0.008; % in um
% Use actual recorded number of datapoints to parse thru data/etc
if any(contains(params, 'N_DATAPOINTS') ~= 0)
   n_datapoints = str2double(values{contains(params, 'N_DATAPOINTS')});
end

% Open motor ID file and read in data 
motor_data_file = fopen(sprintf('%s_motorID.file', sim_name));
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

for i_data = 1 : n_datapoints - 1
    for i_mt = 1:1:n_mts
        motor_IDs = motor_data(:, i_mt, i_data);
        future_IDs = motor_data(:, i_mt, i_data + 1);
        % Do not count motors that are jammed or at the plus end
        jammed_motors = [];
        for i_site = 1 : n_sites
            motor_ID = motor_IDs(i_site);
            if motor_ID == -1
                continue;
            end
            if i_site == 1
                jammed_motors = [jammed_motors motor_ID];
            else
                fwd_ID = motor_IDs(i_site - 1);
                if fwd_ID == -1
                   continue; 
                end
                if fwd_ID ~= motor_ID
                    jammed_motors = [jammed_motors motor_ID];
                end
            end
        end
        %{
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
        %}
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
                delta_time = abs(i_data - start_datapoint);
                run_time = delta_time * time_per_datapoint;
                velocity = (run_length / run_time) * 1000; % convert to nm/s
                % If time bound is above time cutoff, add to data
                %if end_site(1) > endtag_boundary && run_time > 0
                if all(jammed_motors(:) ~= motor_ID)
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
ci_run = paramci(run_dist);
err_run = abs(ci_run(2) - ci_run(1)) / 2;

% matlab's exponential fit always goes to zero; offset it appropriately
min_time = min(lifetimes);
lifetimes = lifetimes - min_time;
% fit distribution of shifted run times
time_dist = fitdist(lifetimes, 'exponential');
% get mean run time
mean_time = time_dist.mu + min_time;
ci_time = paramci(time_dist);
err_time = abs(ci_time(2) - ci_time(1)) / 2;

vel_dist = fitdist(velocities, 'normal');
mean_vel = vel_dist.mu;
ci_vel = paramci(vel_dist);
err_vel = abs(ci_vel(2) - ci_vel(1)) / 2;


mot_stats(1) = mean_run;
mot_stats(2) = err_run;
mot_stats(3) = mean_time;
mot_stats(4) = err_time;
mot_stats(5) = mean_vel;
mot_stats(6) = err_vel;

fprintf('For sim %s (%i runs):\n ', sim_name, n_runs);
disp(mot_stats);

end
