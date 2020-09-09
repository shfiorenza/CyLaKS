clear all;
% Often-changed variables
simName = 'xlink_diffusion_double';
n_sites = [5000, 5000];
max_sites = max(n_sites);
n_mts = length(n_sites);
starting_point = 001;
max_tau = 3.6; % in seconds
% Pseudo-constant variables
delta_t = 0.000025;
n_steps = 60000000;
n_datapoints = 10000;
time_per_datapoint = delta_t * n_steps / n_datapoints;
active_datapoints = n_datapoints - starting_point;
site_size = 0.008; % in um
xlink_cutoff = 4;

fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
xlinkFileStruct = '%s_xlinkID.file';
mtFileStruct = '%s_mt_coord.file';

xlinkFileName = sprintf(fileDirectory, sprintf(xlinkFileStruct, simName));
xlink_data_file = fopen(xlinkFileName);
xlink_raw_data = fread(xlink_data_file, n_mts * max_sites * n_datapoints, '*int');
fclose(xlink_data_file);
xlink_data = reshape(xlink_raw_data, max_sites, n_mts, n_datapoints);

mtFileName = sprintf(fileDirectory, sprintf(mtFileStruct, simName));
mt_data_file = fopen(mtFileName);
mt_raw_data = fread(mt_data_file, n_mts * n_datapoints, '*int');
fclose(mt_data_file);
mt_data = reshape(mt_raw_data, n_mts, n_datapoints);

min_tau = time_per_datapoint;
tau_increment = time_per_datapoint;

taus = min_tau:tau_increment:max_tau;
n_taus = length(taus);

squared_displacements = zeros(n_taus, n_datapoints);
n_entries = zeros(n_taus, 1);

for i_tau = 1:n_taus
    tau = taus(i_tau); % min_tau + (i_tau - 1) * tau_increment;
    tau_step = int32(tau / time_per_datapoint);

    for i_data = (tau_step + starting_point):n_datapoints
        % Data from current time point
        cur_IDs = squeeze(xlink_data(:, :, i_data));
        cur_coords = squeeze(mt_data(:, i_data));
        % Data from previous time point dictated by tau_step
        prev_IDs = squeeze(xlink_data(:, :, i_data - tau_step));
        prev_coords = squeeze(mt_data(:, i_data - tau_step));
        % Scan through ID data at current timestep
        for i_site = 1:max_sites

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
                    i_entry = n_entries(i_tau) + 1;
                    squared_displacements(i_tau, i_entry) = dx^2;
                    n_entries(i_tau) = n_entries(i_tau) + 1;
                    % Otherwise if two MTs, check to see if doubly-bound
                elseif n_mts == 2
                    i_neighb = find(cur_IDs(:, 2) == cur_IDs(i_site, 1), 1);
                    % If i_neighb isn't found, xlink isn't doubly-bound; skip it
                    if isempty(i_neighb)
                        continue;
                    end

                    neighb_coord = i_neighb + cur_coords(2);
                    anchor_coord = (site_coord + neighb_coord) / 2;
                    prev_index = find(prev_IDs(:, 1) == cur_IDs(i_site, 1), 1);
                    % If xlink isn't found in previous timestep, skip
                    if (isempty(prev_index))
                        continue;
                    end

                    prev_site_coord = prev_index + prev_coords(1);
                    prev_i_neighb = find(prev_IDs(:, 2) == prev_IDs(prev_index, 1), 1);
                    % If xlink isn't doubly-bound in prev timestep; skip
                    if (isempty(prev_i_neighb))
                        continue;
                    end

                    prev_neighb_coord = prev_i_neighb + prev_coords(2);
                    prev_anchor_coord = (prev_site_coord + prev_neighb_coord) / 2;
                    dx = double(anchor_coord - prev_anchor_coord) * site_size;
                    i_entry = n_entries(i_tau) + 1;
                    squared_displacements(i_tau, i_entry) = dx^2;
                    n_entries(i_tau) = n_entries(i_tau) + 1;
                else
                    disp('Greater than 2 MTs not implemented yet.');
                    return
                end

            end

        end

    end

end

MSD = zeros(n_taus, 1);
MSD_err = zeros(n_taus, 1);

for i_tau = 1:n_taus
    active_range = 1:n_entries(i_tau);
    MSD(i_tau) = mean(squared_displacements(i_tau, active_range));
    variance = 0.0;

    for i_entry = active_range
        diff_sq = (MSD(i_tau) - squared_displacements(i_tau, i_entry))^2;
        variance = variance + diff_sq / (n_entries(i_tau) - 1);
    end

    MSD_err(i_tau) = sqrt(variance / (n_entries(i_tau)));
end

% Use MATLAB's linear regression to find slope of data
taus = (min_tau:tau_increment:max_tau)';
taus_padded = [ones(length(taus), 1) taus];
fit = taus_padded \ MSD;
% slope = 2*D for MSD plot
D = fit(2) / 2
%}

fig1 = figure();
set(fig1, 'Position', [50, 50, 960, 600])
% Plot data
errorbar(taus, MSD, MSD_err, 'o', 'MarkerSize', 10, 'LineWidth', 2);
hold on
% Plot fit
plot(taus, taus_padded * fit, '--', 'LineWidth', 2)
ylabel('Mean squared displacement (\mum^2)', 'FontSize', 14);
xlabel('Tau (s)', 'FontSize', 14);
title('Doubly-bound crosslinkers', 'FontSize', 16);
legend('Sim data', 'Linear fit', 'location', 'northwest', 'FontSize', 14);
xlim([0.0 (max_tau + tau_increment)]);
ax = gca;
ax.FontSize = 14;
% Plot D_observed obtained via linear regression
dim = [0.142 0.7 .5 .1];
str = sprintf("D_{obs} = %#.3g", D) + " \mum^2s^{-1}";
annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'FontSize', 14);
