%{
clear variables;
fileDirectory = "/home/shane/projects/CyLaKS/run_xlink_diffusion";
simNames = ["xlink_diffusion_single", "xlink_diffusion"];
n_seeds = [11, 4];
%seeds = [0, 1, 2, 3];

n_taus = 10;
tau_increment = 0.25;  %s
time_per_datapoint = 0.01; %s
n_datapoints_max = 100000;
site_size = 0.008; % um
xlink_cutoff = 5;  %sites

%{
exp_tau = [0.15, 0.3, 0.45, 0.60, 0.75, 0.9, 1.1, 1.2, 1.4, 1.5];
exp_msd = [0.007, 0.01, 0.018, 0.016, 0.026, 0.035, 0.038, 0.029, 0.034, 0.030];
exp_err_msd = [0.002, 0.002, 0.004, 0.005, 0.014, 0.013, 0.017, 0.014, 0.016, 0.013];
%}

n_runs = length(simNames);
n_seeds_max = max(n_seeds); %length(seeds);
dx_sq = zeros(n_runs, n_taus, n_datapoints_max*n_seeds_max);
n_entries = zeros(n_runs, n_taus);
MSD = zeros(n_runs, n_taus);
MSD_err = zeros(n_runs, n_taus);

min_tau = time_per_datapoint;
max_tau = (n_taus - 1)*tau_increment + min_tau;
taus = min_tau:tau_increment:max_tau;

for i_run = 1 : n_runs
    for i_seed = 1 : n_seeds(i_run) % n_seeds
        simName = sprintf('%s_%i', simNames(i_run), i_seed - 1)
        % Open log file and parse it into param labels & their values
        log_file = sprintf('%s/%s.log', fileDirectory, simName);
        log = textscan(fileread(log_file), '%s %s', 'Delimiter', '=');
        params = log{1, 1};
        values = log{1, 2};
        % Read in number of MTs
        n_mts = sscanf(values{contains(params, "count ")}, '%g');
        if any(contains(params, "COUNT ") ~= 0)
            n_mts = sscanf(values{contains(params, "COUNT ")}, '%g');
        end
        n_sites = zeros(1, n_mts);
        for i_mt = 1 : n_mts
            string = sprintf("n_sites[%i] ", i_mt - 1);
            n_sites(i_mt) = sscanf(values{contains(params, string)}, '%i');
            if any(contains(params, sprintf("N_SITES[%i] ", i_mt - 1)) ~= 0)
                string = sprintf("N_SITES[%i] ", i_mt - 1);
                n_sites(i_mt) = sscanf(values{contains(params, string)}, '%i');
            end
        end
        % Read in system params
        dt = sscanf(values{contains(params, "dt ")}, '%g');
        steps_per_datapoint = str2double(values{contains(params, "n_steps_per_snapshot ")});
        time_per_datapoint = dt * steps_per_datapoint;
        n_datapoints = str2double(values{contains(params, "n_datapoints ")});
        % Use actual recorded number of datapoints to parse thru data/etc
        if any(contains(params, "N_DATAPOINTS ") ~= 0)
            n_datapoints = str2double(values{contains(params, "N_DATAPOINTS ")});
        end
        max_sites = max(n_sites);
        xlinkFileStruct = '%s_protein_id.file';
        xlinkFileName = sprintf('%s/%s', fileDirectory, sprintf(xlinkFileStruct, simName));
        xlink_data_file = fopen(xlinkFileName);
        xlink_raw_data = fread(xlink_data_file, n_mts * max_sites * n_datapoints, '*int');
        fclose(xlink_data_file);
        xlink_data = reshape(xlink_raw_data, max_sites, n_mts, n_datapoints);
        
        
        mt_data = zeros(n_mts, n_datapoints); %reshape(mt_raw_data, n_mts, n_datapoints);
            
        for i_tau = 1:n_taus
            tau = taus(i_tau); % min_tau + (i_tau - 1) * tau_increment;
            tau_step = int32(tau / time_per_datapoint);
            for i_data = 1 : (n_datapoints - tau_step)
                % Data from current time point
                cur_IDs = squeeze(xlink_data(:, :, i_data));
                cur_coords = squeeze(mt_data(:, i_data));
                % Data from previous time point dictated by tau_step
                prev_IDs = squeeze(xlink_data(:, :, i_data + tau_step));
                prev_coords = squeeze(mt_data(:, i_data + tau_step));
                % Scan through ID data at current timestep
                for i_site = 1 : max_sites
                    if cur_IDs(i_site, 1) ~= -1
                        site_coord = i_site + cur_coords(1);
                        % If one MT, just search over single MT in prev. timestep
                        if n_mts == 1
                            prev_index = find(prev_IDs(:, 1) == cur_IDs(i_site, 1));    
                            if (isempty(prev_index))
                                continue
                            end
                            prev_site_coord = prev_index + prev_coords(1);
                            dx = double(site_coord - prev_site_coord) * site_size;
                            i_entry = n_entries(i_run, i_tau) + 1;
                            dx_sq(i_run, i_tau, i_entry) = dx^2;
                            n_entries(i_run, i_tau) = n_entries(i_run, i_tau) + 1;
                            % Otherwise if two MTs, check to see if doubly-bound
                        elseif n_mts == 2
                            i_neighb = -1;
                            for delta = -xlink_cutoff : xlink_cutoff
                                i_probe = i_site + delta;
                                if i_probe < 1 || i_probe > max_sites
                                    continue
                                end
                                neighb_ID = cur_IDs(i_probe, 2);
                                if neighb_ID == cur_IDs(i_site, 1)
                                    i_neighb = i_probe;
                                end
                            end
                            % If i_neighb isn't found, xlink isn't doubly-bound; skip it
                            if i_neighb == -1
                                continue;
                            end
                            neighb_coord = i_neighb + cur_coords(2);
                            anchor_coord = (site_coord + neighb_coord) / 2;
                            prev_index = find(prev_IDs(:, 1) == cur_IDs(i_site, 1), 1);
                            % If xlink isn't found in previous timestep, skip
                            if (isempty(prev_index))
                                continue;
                            end
                            prev_i_neighb = -1;
                            for delta = -xlink_cutoff : xlink_cutoff
                                i_probe = prev_index + delta;
                                if i_probe < 1 || i_probe > max_sites
                                    continue
                                end
                                neighb_ID = prev_IDs(i_probe, 2);
                                if neighb_ID == prev_IDs(prev_index, 1)
                                    prev_i_neighb = i_probe;
                                end
                            end
                            if prev_i_neighb == -1
                                continue;
                            end
                            prev_site_coord = prev_index + prev_coords(1);
                            prev_neighb_coord = prev_i_neighb + prev_coords(2);
                            prev_anchor_coord = (prev_site_coord + prev_neighb_coord) / 2;
                            dx = double(anchor_coord - prev_anchor_coord) * site_size;
                            %dx = double(i_neighb - prev_i_neighb) * site_size;
                            i_entry = n_entries(i_run, i_tau) + 1;
                            dx_sq(i_run, i_tau, i_entry) = dx^2;
                            n_entries(i_run, i_tau) = n_entries(i_run, i_tau) + 1;
                        else
                            disp('Greater than 2 MTs not implemented yet.');
                            return
                        end
                    end
                end
            end
        end
    end
end
%}

for i_run = 1 : n_runs
    for i_tau = 1 : n_taus
        active_range = 1 : n_entries(i_run, i_tau);
        MSD(i_run, i_tau) = mean(dx_sq(i_run, i_tau, active_range));
        variance = 0.0;
        for i_entry = active_range
            diff_sq = (MSD(i_run, i_tau) - dx_sq(i_run, i_tau, i_entry))^2;
            variance = variance + diff_sq / (n_entries(i_run, i_tau) - 1);
        end
        MSD_err(i_run, i_tau) = sqrt(variance / (n_entries(i_run, i_tau)));
    end
end

% Use MATLAB's linear regression to find slope of data

taus_padded = [ones(length(taus), 1) taus'];
for i_run = 1 : n_runs
    fit = taus_padded \ MSD(i_run, :)';
    % slope = 2*D for MSD plot
    D = fit(2) / 2
    %}
end



fig1 = figure();
set(fig1, 'Position', [50, 50, 720, 540])
hold on
% Plot data
color = [0 0.447 0.741; 0.85, 0.325, 0.098; 0.929, 0.694, 0.125; ...
    0.494, 0.184, 0.556; 0.466, 0.674, 0.188; 0.301, 0.745, 0.933];
D_exp = [0.121, 0.021]; %[0.024, 0.131];
for i_run = 1 : n_runs
    
    plot(taus, 2 * D_exp(i_run) * taus, 'LineWidth', 2, 'Color', [0.6, 0.6, 0.6]);
    %errorbar(taus, MSD, MSD_err, 'o', 'MarkerSize', 10, 'LineWidth', 2);
    sim_data = errorbar(taus, MSD(i_run, :), MSD_err(i_run, :), 'o', 'MarkerSize', 12, 'MarkerEdgeColor', color(i_run, :));
    sim_data.MarkerFaceColor = sim_data.MarkerEdgeColor;
    sim_data.Color = sim_data.MarkerFaceColor;
    % Plot fit
    %plot(taus, taus_padded * fit, 'LineWidth', 2)
    
end
%errorbar(taus, taus * 2 * 0.024,  2 * taus * 2 * 0.003, 'LineWidth', 2);
%errorbar(exp_tau, exp_msd, exp_err_msd);
ylabel('MSD (\mum^2)', 'FontSize', 22);
xlabel('Tau (s)', 'FontSize', 22);
set(gca, 'FontSize', 22);
xlim([-tau_increment (max_tau + tau_increment)]);
ylim([-0.05 0.75]);
xticks([0 1 2]);
yticks([0 0.2 0.4 0.6])
% Force avg_occu plot on top of avg_occu_theory plot
h = get(gca,'Children');
%set(gca,'Children',[h(2) h(1)])
legend([h(3) h(1) h(4)], 'Simulation (1HB)', 'Simulation (2HB)', 'Linear fit', ...
    'location', 'northwest', 'FontSize', 20);
legend('boxoff');

% Plot D_observed obtained via linear regression
%dim = [0.255 0.71 .5 .1];
%str = sprintf("D = %#.3g", D) + " ?m^2s^{-1}";
%annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'FontSize', 18);

%}