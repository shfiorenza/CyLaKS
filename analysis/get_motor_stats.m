function mot_stats = get_motor_stats(base_name, seeds, plot_flag)
% Code block below is for when this is used as a standalone script

%clear variables; 
%base_name = "/home/shane/projects/CyLaKS/test2"; %run_motor_mobility/kif4a_mobility";
%seeds = [0]; %, 1, 2, 3, 4, 5];
%do_plot = true;


do_plot = false;
if nargin == 1 %|| isempty(seeds)
    multirun = false;
    n_seeds = 1;
    sim_name = base_name;
else
    multirun = true;
    n_seeds = length(seeds);
    sim_name = sprintf('%s_%i', base_name, seeds(1));
    %sim_name = base_name;
    if nargin == 3
       do_plot = plot_flag; 
    end
end
%do_plot = true
%}

% Open log file and parse it into param labels & their values
log_file = sprintf('%s.log', sim_name);
log = textscan(fileread(log_file), '%s %s', 'Delimiter', '=');
params = log{1, 1};
values = log{1, 2};
% Read in number of MTs
n_mts = sscanf(values{contains(params, "count ")}, '%g');
if any(contains(params, "COUNT ") ~= 0)
    n_mts = sscanf(values{contains(params, "COUNT ")}, '%g');
end
n_sites = sscanf(values{contains(params, "n_sites[0]")}, '%i');
site_size = 0.008; % in um

runlengths = zeros([(n_seeds*n_mts*n_sites) 1]);
lifetimes = zeros([(n_seeds*n_mts*n_sites) 1]);
velocities = zeros([(n_seeds*n_mts*n_sites) 1]);
n_runs = 0;
for i_seed = 1 : n_seeds
    if multirun
        sim_name = sprintf('%s_%i', base_name, seeds(i_seed));
    else
        sim_name = base_name;
    end
   
    % Open log file and parse it into param labels & their values
    log_file = sprintf('%s.log', sim_name);
    log = textscan(fileread(log_file), '%s %s', 'Delimiter', '=');
    params = log{1, 1};
    values = log{1, 2};
    % Read in system params
    dt = sscanf(values{contains(params, "dt ")}, '%g');
    steps_per_datapoint = str2double(values{contains(params, "n_steps_per_snapshot ")});
    time_per_datapoint = dt * steps_per_datapoint;
    n_datapoints = str2double(values{contains(params, "n_datapoints ")});
    % Use actual recorded number of datapoints to parse thru data/etc
    if any(contains(params, "N_DATAPOINTS ") ~= 0)
        n_datapoints = str2double(values{contains(params, "N_DATAPOINTS ")});
    end
    
    % Open motor ID file and read in data
    motor_data_file = fopen(sprintf('%s_protein_id.file', sim_name));
    raw_motor_data = fread(motor_data_file, n_mts * n_sites * n_datapoints, '*int');
    fclose(motor_data_file);
    motor_data = reshape(raw_motor_data, n_sites, n_mts, n_datapoints);
    
    % have an active list for each MT
    active_motors = zeros([n_mts n_mts*n_sites]);
    n_active = zeros([n_mts 1]);
    
    % motor ID is unique; make following arrays serial w/ ID as index
    starting_site = zeros([10*n_mts*n_sites 1]) - 1;
    starting_datapoint = zeros([10*n_mts*n_sites 1]) - 1;
    
    for i_data = 1 : n_datapoints - 1
        for i_mt = 1:1:n_mts
            motor_IDs = motor_data(:, i_mt, i_data);
            future_IDs = motor_data(:, i_mt, i_data + 1);
            
            % Determine end-tag region; ignore motors that terminate from here
            endtag_boundary = 1;
            for i_site=1:n_sites
                motor_ID = motor_IDs(i_site);
                if motor_ID ~= -1
                   endtag_boundary = i_site + 1;
                else
                    break;
                end
            end
            %}
            % Scan through IDs of bound motors (-1 means no motor on that site)
            for i_site = 1 : n_sites
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
            for i_motor = 1 : 1 : n_active(i_mt)
                i_adj = i_motor - n_deleted;
                motor_ID = active_motors(i_mt, i_adj);
                future_site = find(future_IDs == motor_ID, 1);
                if isempty(future_site)
                    % Calculate distance traveled
                    end_site = find(motor_IDs == motor_ID);
                    start_site = starting_site(motor_ID);
                    delta = abs(end_site(1) - start_site);
                    run_length = delta * site_size * 1000; % convert to nm
                    % Calculate time bound
                    start_datapoint = starting_datapoint(motor_ID);
                    delta_time = abs(i_data - start_datapoint);
                    run_time = delta_time * time_per_datapoint;
                    velocity = run_length / run_time;
                   
                    %{
                    % Determine if we can observe this unbinding event
                    scan_window = 30;
                    observable = true;
                    for i_delta = 1 : 1 : scan_window
                        i_scan_fwd = end_site(1) + i_delta; 
                        if i_scan_fwd < n_sites + 1
                           if future_IDs(i_scan_fwd) ~= -1 && future_IDs(i_scan_fwd) ~= motor_ID
                               observable = false;
                               break;
                           end
                        end 
                        i_scan_bck = end_site(1) - i_delta; 
                        if i_scan_bck > 0
                            if future_IDs(i_scan_bck) ~= -1 && future_IDs(i_scan_bck) ~= motor_ID
                               observable = false; 
                               break;
                            end
                        end
                    end
                    %}
                    
                    % If time bound is above time cutoff, add to data
                    if delta_time > 0 && end_site(1) > endtag_boundary
                        %if all(jammed_motors(:) ~= motor_ID)
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
end
% trim arrays to get rid of un-used containers
runlengths = runlengths(1:n_runs);
lifetimes = lifetimes(1:n_runs);
velocities = velocities(1:n_runs);

mean_run = mean(runlengths);
mean_time = mean(lifetimes);
mean_vel = mean(velocities);
var_run = 0.0;
var_time = 0.0;
var_vel = 0.0;
for i_run = 1 : n_runs
   delta_run = mean_run - runlengths(i_run);
   var_run = var_run + delta_run * delta_run / (n_runs - 1);
   delta_time = mean_time - lifetimes(i_run);
   var_time = var_time + delta_time * delta_time / (n_runs - 1);
   delta_vel = mean_vel - velocities(i_run);
   var_vel = var_vel + delta_vel * delta_vel / (n_runs - 1);
end
err_run = sqrt(var_run / n_runs);
err_time = sqrt(var_time / n_runs);
err_vel = sqrt(var_vel / n_runs);

mot_stats(1) = mean_run;
mot_stats(2) = err_run;
mot_stats(3) = mean_time;
mot_stats(4) = err_time;
mot_stats(5) = mean_vel;
mot_stats(6) = err_vel;

%}
fprintf('For sim %s (%i runs):\n', sim_name, n_runs);
disp(mot_stats);

if do_plot
    % prep figure
    fig1 = figure();
    set(fig1, 'Position', [50, 50, 960, 600]);
    % plot run length histogram
    n_bins = int32(sqrt(n_runs));
    histfit(runlengths, n_bins, 'exponential');
    % Display mean runlength
    dim1 = [0.55 0.65 0.2 0.2];
    str1 = sprintf('Mean processivity: %i +/- %i nm', int32(mean_run), int32(err_run));
    annotation('textbox', dim1, 'String', str1, 'FontSize', 14, 'FitBoxToText', 'on');
    % Display mean lifetime
    dim2 = [0.55 0.58 0.2 0.2];
    str2 = sprintf('Mean life time: %#.2f +/- %#.2f s', mean_time, err_time);
    annotation('textbox', dim2, 'String', str2, 'FontSize', 14, 'FitBoxToText', 'on');
    % Display mean velocity
    dim3 = [0.55 0.51 0.2 0.2];
    str3 = sprintf('Mean velocity: %i +/- %i nm/s', int32(mean_vel), int32(err_vel));
    annotation('textbox', dim3, 'String', str3, 'FontSize', 14, 'FitBoxToText', 'on');
    % Cosmetic stuff
    xlabel('Processivity (nm)', 'FontSize', 20);
    ylabel('Counts', 'FontSize', 20);
    set(gca, 'FontSize', 20);
    
    fig2 = figure();
    set(fig2, 'Position', [75, 75, 960, 600]);
    histfit(lifetimes, n_bins, 'exponential');
    xlabel('Lifetime (s)', 'FontSize', 20);
    ylabel('Counts', 'FontSize', 20);
    set(gca, 'FontSize', 20);
    
    fig3 = figure();
    set(fig3, 'Position', [100, 100, 960, 600]);
    histfit(velocities, n_bins, 'normal');
    xlabel('Velocity (nm/s)', 'FontSize', 20);
    ylabel('Counts', 'FontSize', 20);
    set(gca, 'FontSize', 20);
    
end

end
