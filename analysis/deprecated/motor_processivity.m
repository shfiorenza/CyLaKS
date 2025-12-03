% FIX ME -- very janky implementation of multi PFs at the moment 

clear variables;
seeds = [0, 1, 2, 3, 4, 5];
%seeds = [0];

%file_dir = '../out_final_motMotility';
file_dir = '..';
sim_name_base = 'out_final/shep_0.1nM_10nM_8_1000_0.6kT_3x_5x';
%sim_name_base = 'outputProto5/shep_0.1nM_50nM_8_1000_0.6kT_3x_5x';
%sim_name_base = 'outputProto6/shep_0.1nM_50nM_3_10000_0.6kT_3x_5x';
%sim_name_base = 'final_motility_1nM_1x';
%sim_name_base = 'shep_1nM_100nM_8_1250_0.6kT_3x_5x_0';

output_folder = 'test';


%chosen_SID = xlink_SID;

% Open log file and parse it into param labels & their values
log_file = sprintf('%s/%s', file_dir, sprintf('%s_0.log', sim_name_base));
%log_file = sprintf('%s/%s', file_dir, sprintf('%s.log', sim_name_base));
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
    sim_name = sprintf('%s_%i', sim_name_base, seeds(i_seed))
    %sim_name = sim_name_base;

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
            %{
        % Do not count motors that are jammed or at the plus end
        jammed_motors = [];
        for i_site = 1:n_sites
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
            %}
            % Determine end-tag region; ignore motors that terminate from here
            %{
        endtag_boundary = 1;
        for i_site = 1:n_sites
            motor_ID = motor_IDs(i_site);
            if motor_ID ~= -1
                endtag_boundary = i_site + 1;
            else
                break;
            end

        end
            %}
            endtag_boundary = 400;

            % Scan through IDs of bound motors (-1 means no motor on that site)
            for i_site = 1:1:n_sites
                occupant = occupancy_data(i_site, i_mt, i_data);
                if occupant ~= chosen_SID
                    continue;
                end
                motor_ID = motor_IDs(i_site);
                % Always count motor on first site
                if motor_ID > 0 && i_site == 1
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
                elseif motor_ID > 0 && motor_IDs(i_site - 1) ~= motor_ID
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
                    future_site = find(motor_data(:, j_mt, i_data + 1 ) == motor_ID, 1);
                    if ~isempty(future_site)
                        unbound = false;
                    end
                end
                %future_site = find(future_IDs == motor_ID, 1);
                if unbound
                    % Calculate distance traveled
                    end_site = -1;
                    end_mt = -1;
                    for j_mt = 1 : 1 : n_mts
                        end_site = find(motor_data(:, j_mt, i_data) == motor_ID, 1);
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
                    delta = abs(end_site(1) - start_site);
                    run_length = delta * site_size;
                    % Calculate time bound
                    start_datapoint = starting_datapoint(motor_ID);
                    delta_time = abs(i_data - start_datapoint);
                    run_time = delta_time * time_per_datapoint;
                    velocity = (run_length / run_time) * 1000; % convert to nm/s
                    % If time bound is above time cutoff, add to data
                    if end_site(1) > endtag_boundary && run_time > time_per_datapoint && velocity > 10 && velocity < 2000
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

%}
% trim arrays to get rid of un-used containers
runlengths = runlengths(1:n_runs);
lifetimes = lifetimes(1:n_runs);
velocities = velocities(1:n_runs);

%disp(velocities(find(velocities <= 0)))

avg_run = sum(runlengths) / n_runs
avg_time = sum(lifetimes) / n_runs
avg_vel = sum(velocities) / n_runs
%}
% matlab's exponential fit always goes to zero; offset it appropriately
min_run = min(runlengths);
runlengths = runlengths - min_run;
% fit distribution of shifted run lengths
run_dist = fitdist(runlengths, 'exponential');
mean_run = run_dist.mu + min_run;
conf_inv_run = paramci(run_dist);
sigma_run = abs(conf_inv_run(2) - conf_inv_run(1)) / 2;

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

%min_vel = min(velocities);
%velocities = velocities - (min_vel - 0.0001);
%disp(velocities < 0)

vel_dist = fitdist(velocities, 'lognormal');
mean_vel = exp(vel_dist.mu + (vel_dist.sigma^2)/2);
mean_vel_geo = exp(vel_dist.mu)
sigma_vel_geo = exp(vel_dist.sigma)
%var_vel = (exp(vel_dist.sigma^2) - 1)*exp(2*vel_dist.mu + vel_dist.sigma^2);
%sigma_vel = sqrt(var_vel);
%}
%{
vel_dist = fitdist(velocities, 'normal');
mean_vel = vel_dist.mu;
conf_inv_vel = paramci(vel_dist);
sigma_vel = std(vel_dist)
%}
%{
% prep figure
fig1 = figure();
set(fig1, 'Position', [50, 50, 960, 600]);
set(gcf, 'DefaultAxesFontName', 'Arial');
set(gcf, 'DefaultTextFontName', 'Arial');
% plot run length histogram
n_bins = int32(sqrt(n_runs));
hist = histfit(runlengths, n_bins, 'exponential');
% Display mean runlength
dim1 = [0.55 0.65 0.2 0.2];
str1 = sprintf('Mean run length: %#.1f +/- %#.1g microns', mean_run, sigma_run);
annotation('textbox', dim1, 'String', str1, 'FitBoxToText', 'on');
% Display mean lifetime
dim2 = [0.55 0.6 0.2 0.2];
str2 = sprintf('Mean run time: %#.1f +/- %#.1g seconds', mean_time, sigma_time);
annotation('textbox', dim2, 'String', str2, 'FitBoxToText', 'on');
% Display mean velocity
dim3 = [0.55 0.55 0.2 0.2];
str3 = sprintf('Mean velocity: %#.1f +/- %#.1f nm/s', mean_vel, sigma_vel);
annotation('textbox', dim3, 'String', str3, 'FitBoxToText', 'on');
% Cosmetic stuff
%title(sprintf('Run length histogram for %g micron MT with %i pM Kif4A', int32(n_sites * 0.008), conc));
% sprintf('k on = 0.000242 nM^{-1}s^{-1} | c eff bind = 800,000 nM | k hydro = 100 s^{-1} | k off i = 0.45 s^{-1}')});
xlabel('Run length (um)');
ylabel('Counts');

fig2 = figure();
set(fig2, 'Position', [75, 75, 960, 600]);
set(gcf, 'DefaultAxesFontName', 'Arial');
set(gcf, 'DefaultTextFontName', 'Arial');
histfit(lifetimes, n_bins, 'exponential');
xlabel('Lifetime (s)');
ylabel('Counts');
fig3 = figure();
set(fig3, 'Position', [100, 100, 720, 600]);
set(gcf, 'DefaultAxesFontName', 'Arial');
set(gcf, 'DefaultTextFontName', 'Arial');
hold on
histfit(velocities, n_bins, 'normal');
xlabel('Velocity (nm/s)');
ylabel('Counts');
%}
% saveas(fig1, sprintf('%s/%s_proc.png', output_folder, sim_name), 'png');
% saveas(fig1, sprintf('%s/%s_proc.svg', output_folder, sim_name), 'svg');
% saveas(fig2, sprintf('%s/%s_time.png', output_folder, sim_name), 'png');
% saveas(fig2, sprintf('%s/%s_time.svg', output_folder, sim_name), 'svg');
% saveas(fig3, sprintf('%s/%s_vel.png', output_folder, sim_name), 'png');
% saveas(fig3, sprintf('%s/%s_vel.svg', output_folder, sim_name), 'svg');



n_bins = 45;
fig3 = figure();
set(fig3, 'Position', [100, 100, 720, 600]);
hold on
h3 = histfit(velocities, n_bins, 'lognormal'); %, 'FaceColor', 'k', 'EdgeColor', [0 0 0], 'LineWidth', 2);
%h3 = histfit(velocities);
h3(1).FaceColor = [0 0 0];
h3(1).EdgeColor = 'none';
h3(1).BarWidth = 0.5;
h3(2).LineWidth = 3;
%title(sprintf('Velocity histogram for %g micron MT with %i pM Kif4A', int32(n_sites * 0.008), conc));
xlabel('Velocity (nm/s)');
ylabel('Frequency');
xlim([0 400]);
%xlim([0 1000]);
dim3 = [0.525 0.65 0.2 0.2];
str3 = sprintf('Geo Mean = %#.1f nm/s', mean_vel_geo);
%str3 = sprintf('Mean = %#.1f nm/s', mean_vel);
annotation('textbox', dim3, 'String', str3, 'FitBoxToText', 'on', 'EdgeColor','none', 'FontSize', 24);
dim4 = [0.525 0.6 0.2 0.2];
str4 = sprintf('Geo SD = %#.1f', sigma_vel_geo);
%str4 = sprintf('SD = %#.1f', sigma_vel);
annotation('textbox', dim4, 'String', str4, 'FitBoxToText', 'on','EdgeColor','none', 'FontSize', 24);
dim4 = [0.525 0.55 0.2 0.2]; % was 0.45
str4 = sprintf('N = %i', length(velocities));
annotation('textbox', dim4, 'String', str4, 'FitBoxToText', 'on','EdgeColor','none', 'FontSize', 24);
set(gca,'box','off')
set(gca, 'FontSize', 28);
set(gca,'TickDir','out');
set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);
%}