
clear all;
applied_forces = [-12, -10, -8, -6, -4, -2, 0];
n_sites = 5000;
n_steps = 12500000;
n_mts = 1;
delta_t = 0.0002; %0.000025;
n_datapoints = 10000;
time_per_datapoint = delta_t * n_steps / n_datapoints;
starting_point = 1;
active_datapoints = n_datapoints - starting_point;
time_cutoff = time_per_datapoint; %0.3;    %in seconds
site_size = 8; % in nm

fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
motorFileStruct = '%s_motorID.file';

exp_forces = [0, -1.0 : -0.5 : -6.0];
exp_runlengths = [940, 770, 480, 500, 270, 330, 250, 150, 240, 120, 70, 50];
exp_err_runlengths = [130, 80, 40, 60, 20, 30, 20, 10, 10, 9, 9, 10];
exp_velocities = [750, 690, 720, 640, 630, 540, 440, 390, 250, 230, 140, 120];
exp_err_velocities = [6, 5, 5, 6, 10, 6, 5, 6, 5, 6, 5, 5];

n_runs = length(applied_forces);
avg_runlengths = zeros(n_runs, 1);
err_runlengths = zeros(n_runs, 1);
avg_lifetimes = zeros(n_runs, 1);
err_lifetimes = zeros(n_runs, 1);
avg_velocities = zeros(n_runs, 1);
err_velocities = zeros(n_runs, 1);
for i_run = 1 : n_runs 
    simName = sprintf("processivity_%ipN", abs(applied_forces(i_run)));
    motorFileName = sprintf(fileDirectory, sprintf(motorFileStruct, simName));
    motor_data_file = fopen(motorFileName);
    raw_motor_data = fread(motor_data_file, [n_mts * n_sites * n_datapoints], '*int');
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
                    % Store runlength, lifetime, & velocity data
                    n_runs = n_runs + 1;
                    runlengths(n_runs) = run_length;
                    lifetimes(n_runs) = run_time;
                    velocities(n_runs) = (run_length / run_time);
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
    conf_inv_run = paramci(run_dist);
    sigma_run = abs(conf_inv_run(2) - conf_inv_run(1)) / 2;
    
    avg_runlengths(i_run) = mean_run;
    err_runlengths(i_run) = sigma_run;
    
    %{
    % matlab's exponential fit always goes to zero; offset it appropriately
    min_time = min(lifetimes);
    lifetimes = lifetimes - min_time;
    % fit distribution of shifted run times
    time_dist = fitdist(lifetimes, 'exponential');
    % get mean run time
    mean_time = time_dist.mu + min_time;
    % get confidence interval of mean run time
    conf_inv_time = paramci(time_dist);
    % calculate sigma
    sigma_time = abs(conf_inv_time(2) - conf_inv_time(1)) / 2;
    %}
    
    vel_dist = fitdist(velocities, 'normal');
    mean_vel = vel_dist.mu;
    conf_inv_vel = paramci(vel_dist);
    sigma_vel = round(abs(conf_inv_vel(2) - conf_inv_vel(1)) / 2, 1, 'significant');
    
    avg_velocities(i_run) = mean_vel;
    err_velocities(i_run) = sigma_vel;
    
end
%}

% Run length plot %
fig1 = figure();
set(fig1, 'Position', [50, 50, 960, 600]);
hold all
% Plot sim data
sim_run = errorbar(applied_forces, avg_runlengths, err_runlengths, 'o', ...
    'MarkerSize', 8, 'LineWidth', 2);
sim_run.MarkerFaceColor = sim_run.MarkerEdgeColor;
% Plot experimental data
exp_run = errorbar(exp_forces, exp_runlengths, exp_err_runlengths, 'sq', ...
    'MarkerSize', 8, 'LineWidth', 2);
exp_run.MarkerFaceColor = exp_run.MarkerEdgeColor;
% Plot experimental fit 
L_0 = 1120;         % in nm
sigma_off = 2.0;    % in nm
kbT = 4.114;        % in pN * nm
exp_run_fit = fplot(@(f) L_0*exp(-abs(f) * sigma_off/kbT), [-6 0], 'LineWidth', 2);
exp_run_fit.Color = exp_run.Color;
% Bring sim data to front
uistack(sim_run, 'top');
ax = gca;
ax.FontSize = 12; 
% Label axes, legend, etc. 
xlabel('Applied force (pN)', 'FontSize', 14);
ylabel('Run length (nm)', 'FontSize', 14);
%xlim([-6.5 0.5]);
legend({'Experimental data', 'Experimental fit', 'Simulation data'}, ... 
    'location', 'northwest', 'FontSize', 12);


fig2 = figure();
set(fig2, 'Position', [75, 75, 960, 600]);
hold all
sim_vel = errorbar(applied_forces, avg_velocities, err_velocities, 'o', ...
    'MarkerSize', 8, 'LineWidth', 2);
sim_vel.MarkerFaceColor = sim_vel.MarkerEdgeColor;
exp_vel = errorbar(exp_forces, exp_velocities, exp_err_velocities, 'sq', ...
    'MarkerSize', 8, 'LineWidth', 2);
exp_vel.MarkerFaceColor = exp_vel.MarkerEdgeColor;
d_step = 8.2;         % in nm
F_i = 26;           % in pN
k0_1 = 4900;         % in 1/s
sigma_1 = 4.6;      % in nm
k0_2 = 95;          % in 1/s
k0_3 = 260;          % in 1/s
sigma_3 = 0.35;     % in nm
k_1 = @(f) k0_1 * exp(f * sigma_1 / kbT);
k_2 = k0_2; 
k_3 = @(f) k0_3 * exp((f + F_i) * sigma_3 / kbT);
exp_vel_fit = fplot(@(f) d_step*k_1(f)*k_2*k_3(f) / (k_1(f)*k_2 + k_3(f)*(k_1(f) + k_2)), ...
    [-6 0], 'LineWidth', 2);
exp_vel_fit.Color = exp_vel.Color;
uistack(sim_vel, 'top');
ax = gca;
ax.FontSize = 12; 
xlabel('Applied force (pN)', 'FontSize', 14);
ylabel('Velocity (nm/s)', 'FontSize', 14);
%xlim([-6.5 0.5]);
legend({'Experimental data', 'Experimental fit', 'Simulation data'}, ... 
    'location', 'northwest', 'FontSize', 12);
