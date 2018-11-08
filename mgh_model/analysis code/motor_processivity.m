
clear all;
% Often-changed variables
n_sites = 1250;
simName = 'mot500kb';
% Pseudo-constant variables
n_mts = 1;
delta_t = 0.00005;
n_steps = 20000000;
n_datapoints = 100000;
time_per_datapoint = delta_t * n_steps / n_datapoints;
starting_point = 1;
active_datapoints = n_datapoints - starting_point;
motor_speed = 0.65;     %in um/s
smallest_time = 0.3;    %in seconds
site_size = 0.008; % in um
smallest_run = motor_speed * smallest_time;  % in um

fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
motorFileStruct = '%s_motorID.file';

conc_eff = [800,1600,8000,16000,80000,160000,800000,1600000, ...
            8000000,16000000,80000000,160000000];

means = zeros([length(conc_eff) 1]);

%for i_conc=1:1:length(conc_eff)
    
    %simName = sprintf('MotorProScan_%d', conc_eff(i_conc));
    
    motorFileName = sprintf(fileDirectory, sprintf(motorFileStruct, simName));
    
    motor_data_file = fopen(motorFileName);
    raw_motor_data = fread(motor_data_file, [n_mts * n_sites * n_datapoints], '*int');
    fclose(motor_data_file);
    motor_data = reshape(raw_motor_data, n_sites, n_mts, n_datapoints);
    
    run_lengths = zeros([(n_mts*n_sites*n_datapoints) 1]);
    n_runs = 0;
    active_motors = zeros([n_mts*n_sites 1]);
    n_active = 0;
    starting_site = zeros([n_mts n_mts*n_sites]) - 1;
    
    for i_data = starting_point:1:n_datapoints - 1
        for i_mt = 1:1:n_mts
            motor_IDs = motor_data(:, i_mt, i_data);
            % Scan through IDs of bound motors (-1 means no motor on that site)
            for i_site = 1:1:n_sites - 1
                motor_ID = motor_IDs(i_site);
                % If a motor is found, check if it's already been accounted for
                if motor_ID > 0 ...
                        && motor_IDs(i_site + 1) ~= motor_ID
                    % Record the motor's starting site if this is the first time
                    % seeing it (-1 means it was not seen last datapoint)
                    if starting_site(i_mt, motor_ID) == -1
                        starting_site(i_mt, motor_ID) = i_site;
                        n_active = n_active + 1;
                        active_motors(n_active) = motor_ID;
                    end
                end
            end
            % Check one datapoint into the future to see if any motors unbound
            future_IDs = motor_data(:, i_mt, i_data + 1);
            n_deleted = 0;
            for i_motor = 1:1:n_active
                i_adj = i_motor - n_deleted;
                motor_ID = active_motors(i_adj);
                future_site = find(future_IDs == motor_ID);
                if isempty(future_site)
                    end_site = find(motor_IDs == motor_ID);
                    start_site = starting_site(i_mt, motor_ID);
                    delta = abs(end_site(1) - start_site);
                    run_length = delta * site_size;
                    if run_length >= smallest_run
                        n_runs = n_runs + 1;
                        run_lengths(n_runs) = run_length;
                    end
                    starting_site(i_mt, motor_ID) = -1;
                    active_motors(i_adj) = -1;
                    active_motors(i_adj:n_active) = circshift(active_motors(i_adj:n_active), -1);
                    n_active = n_active - 1;
                    n_deleted = n_deleted + 1;
                end
            end
        end
    end
    
    run_lengths = run_lengths(1:n_runs);
    n_bins = int32(sqrt(length(run_lengths)));
    % matlab's exponential fit always goes to zero...offset it appropriately
    min_val = min(run_lengths);
    run_lengths = run_lengths - min_val;
    
    fig1 = figure();
    set(fig1, 'Position', [50, 50, 960, 600])
    
        % run histfit with exponential fit
    hist = histfit(run_lengths, n_bins, 'exponential');
    
    title('10 micron long MT; c\_motor = 0.5 nM; c\_eff\_bind = 700000 nM'); 
    xlabel('Run length (um)');
    ylabel('Counts');
   
    % get x/y data of fit
    x=get(hist(2),'XData');
    y=get(hist(2),'YData');
    % normalize fit
    max_height = max(y);
    y = y / max_height;
    % find value and index of entry in y that is closest to mean
    [val, index] = min(abs(y - exp(-1)));
    % get un-normalized value of mean by scaling
    mean = double(index) / double(length(y)) * x(length(x)) + min_val;
    % add min value back to x
    x = x + min_val;
  
    dim = [0.55 0.75 .1 .1];
    str = sprintf('Mean run length: %#.3g microns', mean);
    annotation('textbox',dim,'String',str,'FitBoxToText','on');  
    
    %means(i_conc) = mean;
%end

%plot(conc_eff, means)

fig1 = figure();
set(fig1, 'Position', [50, 50, 960, 600])
title('10 micron long MT; c\_motor = 0.5 nM'); 
xlabel('Effective binding concentration of 2nd head (nM)');
ylabel('Mean processivity (um)');
grid on
grid minor
%semilogx(conc_eff, means, 'LineWidth', 2)
