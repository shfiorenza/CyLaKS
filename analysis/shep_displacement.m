% FIX ME -- very janky implementation of multi PFs at the moment 

clear variables;
%seeds = [0, 1, 2, 3, 4, 5];
seeds = [0];

file_dir = '..';

%sim_name_base = 'out_final_xlinkOnly/shep_0.1nM_0.0nM_8_1000_0.6kT_3x_5x_0';
%sim_name_base = 'out_final/shep_1nM_100nM_8_1000_0.6kT_3x_5x_0';
%sim_name_base = 'out_final_xlinkOnly/shep_0.1nM_0.0nM_8_1000_0.6kT_3x_5x_0';
sim_name_base = 'out_final/shep_0.1nM_100nM_8_1000_0.6kT_3x_5x_0';
%sim_name_base = 'out_final/shep_0.1nM_10nM_8_1000_0.6kT_3x_5x_0'
sim_name_base = '/out_final_xlinkLifetime/shep_0.01nM_100nM_8_1000_0.6kT_3x_5x_0_xlink_0.1x'

xlink_SID = 1;
motor_SID = 2;
chosen_SID = motor_SID; %xlink_SID;

% Open log file and parse it into param labels & their values
%log_file = sprintf('%s/%s', file_dir, sprintf('%s_0.log', sim_name_base));
log_file = sprintf('%s/%s', file_dir, sprintf('%s.log', sim_name_base));
log = textscan(fileread(log_file), '%s %s', 'Delimiter', '=');
params = log{1, 1};
values = log{1, 2};
% Read in number of MTs
n_mts = sscanf(values{contains(params, "count ")}, '%g');
if any(contains(params, "n_subfilaments") ~= 0)
    n_sub = sscanf(values{contains(params, "n_subfilaments ")}, '%g');
    if n_sub > n_mts
       n_mts = n_sub;
    end
end
if any(contains(params, "COUNT ") ~= 0)
    n_mts = sscanf(values{contains(params, "COUNT ")}, '%g');
end
mt_lengths = zeros(1, n_mts);
for i_mt = 1 : n_mts
    string = sprintf("n_sites[%i] ", i_mt - 1);
    mt_lengths(i_mt) = sscanf(values{contains(params, string)}, '%i');
    if any(contains(params, sprintf("N_SITES[%i] ", i_mt - 1)) ~= 0)
        string = sprintf("N_SITES[%i] ", i_mt - 1);
        mt_lengths(i_mt) = sscanf(values{contains(params, string)}, '%i');
    end
end
n_sites = max(mt_lengths);
% Read in system params
dt = sscanf(values{contains(params, "dt ")}, '%g');
steps_per_datapoint = str2double(values{contains(params, "n_steps_per_snapshot ")});
time_per_datapoint = dt * steps_per_datapoint;
n_datapoints = str2double(values{contains(params, "n_datapoints ")});
% Use actual recorded number of datapoints to parse thru data/etc
if any(contains(params, "N_DATAPOINTS ") ~= 0)
    n_datapoints = str2double(values{contains(params, "N_DATAPOINTS ")});
end
n_dims = 2;
site_size = 0.0082; % in um
n_seeds = length(seeds);

% motor ID is unique; make following arrays serial w/ ID as index
runlengths = zeros([(100 * n_seeds * n_mts * n_sites) 1]);
lifetimes = zeros([(100 * n_seeds  * n_mts * n_sites) 1]);
velocities = zeros([(100 * n_seeds  * n_mts * n_sites) 1]);
n_runs = 0;

for i_seed = 1:n_seeds
    %sim_name = sprintf('%s_%i', sim_name_base, seeds(i_seed))
    sim_name = sim_name_base;

    motorFileStruct = '%s_protein_id.file';
    motorFileName = sprintf("%s/%s", file_dir, sprintf(motorFileStruct, sim_name));
    motor_data_file = fopen(motorFileName);
    raw_motor_data = fread(motor_data_file, n_mts * n_sites * n_datapoints, '*int');
    fclose(motor_data_file);
    motor_data = reshape(raw_motor_data, n_sites, n_mts, n_datapoints);

    occuFileStruct = '%s_occupancy.file';
    occuFileName = sprintf("%s/%s", file_dir, sprintf(occuFileStruct, sim_name));
    occu_data_file = fopen(occuFileName);
    raw_occu_data = fread(occu_data_file, n_mts * n_sites * n_datapoints, '*int');
    fclose(occu_data_file);
    occupancy_data = reshape(raw_occu_data, n_sites, n_mts, n_datapoints);

    % have an active list for each MT
    active_motors = zeros([n_mts n_mts * n_sites * 10]);
    n_active = zeros([n_mts 1]);

    starting_site = zeros([100 * n_mts * n_sites 1]) - 1;
    starting_mt = zeros([100 * n_mts * n_sites 1]) - 1;
    starting_datapoint = zeros([100 * n_mts * n_sites 1]) - 1;
    for i_data = 1 : n_datapoints - 1
        for i_mt = 1:1:n_mts
            motor_IDs = motor_data(:, i_mt, i_data);
            %future_IDs = motor_data(:, i_mt, i_data + 1);
            endtag_boundary = 125;

            % Scan through IDs of bound motors (-1 means no motor on that site)
            for i_site = 1:1:n_sites
                occupant = occupancy_data(i_site, i_mt, i_data);
                if occupant ~= chosen_SID
                    continue;
                end
                motor_ID = motor_IDs(i_site);
                % Always count motor on first site
                if motor_ID > 0 && i_site == n_sites
                    % Record the motor's starting site if this is the first time
                    % seeing it (-1 means it was not seen last datapoint)
                    if starting_site(motor_ID) == -1
                        starting_site(motor_ID) = i_site;
                        starting_datapoint(motor_ID) = i_data;
                        starting_mt(motor_ID) =  i_mt;
                        n_active(i_mt) = n_active(i_mt) + 1;
                        active_motors(i_mt, n_active(i_mt)) = motor_ID;
                    end
                    % Otherwise if a motor is found, only count first head
                elseif motor_ID > 0 && motor_IDs(i_site + 1) ~= motor_ID
                    % Record the motor's starting site if this is the first time
                    % seeing it (-1 means it was not seen last datapoint)
                    if starting_site(motor_ID) == -1
                        starting_site(motor_ID) = i_site;
                        starting_datapoint(motor_ID) = i_data;
                        starting_mt(motor_ID) =  i_mt;
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
                unbound = true;
                for j_mt = 1 : 1 : n_mts
                    future_site = find(motor_data(:, j_mt, i_data + 1 ) == motor_ID);
                    if ~isempty(future_site)
                        if chosen_SID == motor_SID
                            if j_mt == starting_mt(motor_ID) && future_site(1) < starting_site(motor_ID)
                                unbound = false;
                            end
                        else
                            unbound = false;
                        end
                    end
                end
                %future_site = find(future_IDs == motor_ID, 1);
                if unbound
                    % Calculate distance traveled
                    end_site = -1;
                    end_mt = -1;
                    for j_mt = 1 : 1 : n_mts
                        end_site = find(motor_data(:, j_mt, i_data) == motor_ID);
                        if ~isempty(end_site)
                            end_mt = j_mt;
                            break;
                        end
                    end
                    if isempty(end_site)
                        disp('Error');
                        return;
                    end
                    start_site = starting_site(motor_ID);
                    delta = end_site(1) - start_site; % we want actual displacement, not MSD
                    run_length = delta * site_size;
                    % Calculate time bound
                    start_datapoint = starting_datapoint(motor_ID);
                    delta_time = i_data - start_datapoint;
                    run_time = delta_time * time_per_datapoint;
                    velocity = (abs(run_length) / run_time) * 1000; % make vel pos. only and convert to nm/s
                    % If time bound is above time cutoff, add to data
                    if end_site(1) > endtag_boundary && run_time > time_per_datapoint && velocity > 10 && velocity < 50000
                        %if run_length > 0
                        %    fprintf("%g, %g\n", run_time, velocity);
                        %end
                        n_runs = n_runs + 1;
                        runlengths(n_runs) = run_length;
                        lifetimes(n_runs) = run_time;
                        velocities(n_runs) = velocity;
                    end
                    starting_site(motor_ID) = -1;
                    starting_datapoint(motor_ID) = -1;
                    starting_mt(motor_ID) = -1;
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
%}

%disp(velocities(find(velocities <= 0)))

avg_run = sum(runlengths) / n_runs;
err_run = std(runlengths) / sqrt(n_runs);
avg_time = sum(lifetimes) / n_runs;
err_time = std(lifetimes) / sqrt(n_runs);
avg_vel = sum(velocities) / n_runs;
err_vel = std(velocities) / sqrt(n_runs);
fprintf("Run: %.2f ± %.2f | Time: %.2f ± %.2f | Vel: %.2f ± %.2f\n", avg_run, err_run, avg_time, err_time, avg_vel, err_run)
%}

% runlengths - normal
run_dist = fitdist(runlengths, 'normal');
mean_run = run_dist.mean;
sigma_run = run_dist.std;
if(abs(avg_run - mean_run) > 1e-6)
    disp("Error in runlength calculations")
    return
end
conf_inv_run = paramci(run_dist);
stderr_run = (conf_inv_run(2) - conf_inv_run(1)) / 4; % get std_err
% lifetimes - exponential
min_time = min(lifetimes);
lifetimes_adj = lifetimes - min_time; % exp. function always goes to zero; need to adjust
time_dist = fitdist(lifetimes_adj, 'exponential');
mean_time = time_dist.mean + min_time;
sigma_time = time_dist.std;
if(abs(avg_time - mean_time) > 1e-6)
    disp("Error in lifetime calculations")
    return
end
conf_inv_time = paramci(time_dist);
stderr_time = (conf_inv_time(2) - conf_inv_time(1)) / 4; % get std_err
% velocities - lognormal
vel_dist = fitdist(velocities, 'lognormal');
mean_vel = vel_dist.mean;
sigma_vel = vel_dist.std;
conf_inv_vel = paramci(vel_dist);
stderr_vel = (conf_inv_vel(2) - conf_inv_vel(1)) / 4; % get std_err
mean_vel_geo = exp(vel_dist.mu);
sigma_vel_geo = exp(vel_dist.sigma);

% plot figures
n_bins = int32(sqrt(n_runs));
% runlengths
runs = figure('Position', [50, 50, 720, 600]);
h = histfit(runlengths, n_bins, 'normal');
% Display fit statistics
dim1 = [0.625 0.7 0.2 0.2];
str1 = sprintf('Mean = %#.2f ± %#.2f um', mean_run, stderr_run);
dim2 = [0.625 0.65 0.2 0.2];
str2 = sprintf('SD = %#.2f um', sigma_run);
dim3 = [0.625 0.6 0.2 0.2];
str3 = sprintf('N = %i', n_runs);
annotation('textbox', dim1, 'String', str1, 'FitBoxToText', 'on', 'EdgeColor','none', 'FontSize', 14);
annotation('textbox', dim2, 'String', str2, 'FitBoxToText', 'on', 'EdgeColor','none', 'FontSize', 14);
annotation('textbox', dim3, 'String', str3, 'FitBoxToText', 'on', 'EdgeColor','none', 'FontSize', 14);
% Cosmetic stuff
title(sprintf('%s', sim_name),'Interpreter', 'none');
xlabel('Run length (um)');
ylabel('Counts');
set(gca,'box','off')
set(gca, 'FontSize', 14);
set(gca,'TickDir','out');
set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);
h(1).FaceColor = [0 0 0];
h(1).EdgeColor = 'none';
h(1).BarWidth = 0.5;
h(2).LineWidth = 3;

%lifetimes 
lifetime = figure('Position', [75, 75, 720, 600]);
h = histfit(lifetimes, n_bins, 'exponential');
% Display fit statistics
dim1 = [0.625 0.7 0.2 0.2];
str1 = sprintf('Mean = %#.2f ± %#.2f s', mean_time, stderr_time);
dim2 = [0.625 0.65 0.2 0.2];
str2 = sprintf('SD = %#.2f s', sigma_time);
dim3 = [0.625 0.6 0.2 0.2];
str3 = sprintf('N = %i', n_runs);
annotation('textbox', dim1, 'String', str1, 'FitBoxToText', 'on', 'EdgeColor','none', 'FontSize', 14);
annotation('textbox', dim2, 'String', str2, 'FitBoxToText', 'on', 'EdgeColor','none', 'FontSize', 14);
annotation('textbox', dim3, 'String', str3, 'FitBoxToText', 'on', 'EdgeColor','none', 'FontSize', 14);
% Cosmetic stuff
title(sprintf('%s', sim_name),'Interpreter', 'none');
xlabel('Lifetime (s)');
ylabel('Counts');
set(gca,'box','off')
set(gca, 'FontSize', 14);
set(gca,'TickDir','out');
set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);
h(1).FaceColor = [0 0 0];
h(1).EdgeColor = 'none';
h(1).BarWidth = 0.5;
h(2).LineWidth = 3;

% velocities
vel = figure('Position', [100, 100, 720, 600]);
h = histfit(velocities, n_bins, 'lognormal');
% Display fit statistics
dim1 = [0.625 0.7 0.2 0.2];
str1 = sprintf('Mean = %#.2f ± %#.2f nm/s', mean_vel, stderr_vel);
dim2 = [0.625 0.65 0.2 0.2];
str2 = sprintf('SD = %#.2f nm/s', sigma_vel);
dim3 = [0.625 0.6 0.2 0.2];
str3 = sprintf('Geo Mean = %#.2f nm/s', mean_vel_geo);
dim4 = [0.625 0.55 0.2 0.2];
str4 = sprintf('Geo SD = %#.2f nm/s', sigma_vel_geo);
dim5 = [0.625 0.5 0.2 0.2];
str5 = sprintf('N = %i', n_runs);
annotation('textbox', dim1, 'String', str1, 'FitBoxToText', 'on', 'EdgeColor','none', 'FontSize', 14);
annotation('textbox', dim2, 'String', str2, 'FitBoxToText', 'on', 'EdgeColor','none', 'FontSize', 14);
annotation('textbox', dim3, 'String', str3, 'FitBoxToText', 'on', 'EdgeColor','none', 'FontSize', 14);
annotation('textbox', dim4, 'String', str4, 'FitBoxToText', 'on', 'EdgeColor','none', 'FontSize', 14);
annotation('textbox', dim5, 'String', str5, 'FitBoxToText', 'on', 'EdgeColor','none', 'FontSize', 14);
% Cosmetic stuff
title(sprintf('%s', sim_name),'Interpreter', 'none');
xlabel('Velocity (nm/s)');
ylabel('Counts');
xlim([0 max(velocities)]);
set(gca,'box','off')
set(gca, 'FontSize', 14);
set(gca,'TickDir','out');
set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);
h(1).FaceColor = [0 0 0];
h(1).EdgeColor = 'none';
h(1).BarWidth = 0.5;
h(2).LineWidth = 3;

%{
saveas(fig1, sprintf('%s_proc.png', sim_name), 'png');
saveas(fig1, sprintf('%s_proc.svg', sim_name), 'svg');
saveas(fig2, sprintf('%s_time.png', sim_name), 'png');
saveas(fig2, sprintf('%s_time.svg', sim_name), 'svg');
saveas(fig3, sprintf('%s_vel.png', sim_name), 'png');
saveas(fig3, sprintf('%s_vel.svg', sim_name), 'svg');
%}
