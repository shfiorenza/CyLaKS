
clear all;
% Often-changed variables
i_conc = 6;
kif4a_concs = [20, 50, 80, 120, 220, 420];
mt_lengths = [75000, 50000, 40000, 25000, 8500, 3250];
conc = kif4a_concs(i_conc);
%n_sites = mt_lengths(i_conc);
n_sites = 5000;
%simName = sprintf("kif4a_coop_optimization_lifetimeOnly_10.1_%i", conc);
%simName = sprintf("lattice_coop_%i", conc);
simName = "test";
% Pseudo-constant variables
n_mts = 1;
delta_t = 0.00002; %0.000025;
n_steps = 12500000;
n_datapoints = 10000;
time_per_datapoint = delta_t * n_steps / n_datapoints;
starting_point = 1;
active_datapoints = n_datapoints - starting_point;
time_cutoff = time_per_datapoint; %0.3;    %in seconds
site_size = 0.008; % in um

fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
motorFileStruct = '%s_motorID.file';

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
        jammed_region = 0;
        jam_start = -1;
        n_jammed = 0;
        %endtag_boundary = 2;
        % Determine end-tag region; ignore motors that terminate here
        for i_site = 1 : n_sites
            motor_ID = motor_IDs(i_site);
            if motor_ID ~= -1
                if jam_start == -1
                    jam_start = i_site;
                end
                n_jammed = n_jammed + 1;
                %endtag_boundary = i_site + 1;
            else
                if n_jammed > 5
                    jammed_region = [jammed_region jam_start:i_site];
                end
                n_jammed = 0;
                jam_start = -1;  
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
                if true %isempty(find(jammed_region == end_site, 1))
                    n_runs = n_runs + 1;
                    runlengths(n_runs) = run_length;
                    lifetimes(n_runs) = run_time;
                    velocities(n_runs) = (run_length / run_time) * 1000; % convert to nm/s  
                    %{
                    if run_time > 100
                       fprintf("Motor %i had lifetime %g\n", motor_ID, run_time); 
                       fprintf("End site: %i, endtag_boundary: %i\n", end_site(1), endtag_boundary);
                    end
                    if velocities(n_runs) < 40
                       fprintf("Motor %i had vel %g\n", motor_ID, velocities(n_runs)); 
                       fprintf("End site: %i, endtag_boundary: %i\n", end_site(1), endtag_boundary);
                    end
                    %}
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

vel_dist = fitdist(velocities, 'normal');
mean_vel = vel_dist.mu;
conf_inv_vel = paramci(vel_dist);
sigma_vel = round(abs(conf_inv_vel(2) - conf_inv_vel(1)) / 2, 1, 'significant');

% prep figure
fig1 = figure();
set(fig1, 'Position', [50, 50, 960, 600]);
% plot run length histogram
n_bins = int32(sqrt(n_runs));
hist = histfit(runlengths, n_bins, 'exponential');
% Display mean runlength
dim1 = [0.55 0.65 0.2 0.2];
str1 = sprintf('Mean run length: %#.1f +/- %#.1g microns', mean_run, sigma_run);
annotation('textbox',dim1,'String',str1,'FitBoxToText','on');
% Display mean lifetime
dim2 = [0.55 0.6 0.2 0.2];
str2 = sprintf('Mean run time: %#.1f +/- %#.1g seconds', mean_time, sigma_time);
annotation('textbox',dim2,'String',str2,'FitBoxToText','on');
% Display mean velocity
dim3 = [0.55 0.55 0.2 0.2];
str3 = sprintf('Mean velocity: %#.1f +/- %#d nm/s', mean_vel, sigma_vel);
annotation('textbox',dim3,'String',str3,'FitBoxToText','on');
% Cosmetic stuff
%title(sprintf('Run length histogram for %g micron MT with %i pM Kif4A', int32(n_sites * 0.008), conc));
   % sprintf('k on = 0.000242 nM^{-1}s^{-1} | c eff bind = 800,000 nM | k hydro = 100 s^{-1} | k off i = 0.45 s^{-1}')});
xlabel('Run length (um)');
ylabel('Counts');

%{
fig2 = figure();
set(fig2, 'Position', [75, 75, 960, 600]);
histfit(lifetimes, n_bins, 'exponential');
%}

fig3 = figure();
set(fig3, 'Position', [100, 100, 960, 600]);
histfit(velocities, n_bins, 'normal');
%title(sprintf('Velocity histogram for %g micron MT with %i pM Kif4A', int32(n_sites * 0.008), conc));
xlabel('Velocity (nm/s)');
ylabel('Counts');