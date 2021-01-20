%
% Often-changed variables
fileDirectory = '/home/shane/projects/CyLaKS/%s';
simName = 'testin2';
n_taus = 10;
i_tau = 0.01;
tau_increment = i_tau;

% Open log file and parse it into param labels & their values
log_file = sprintf(fileDirectory, sprintf('%s.log', simName));
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
site_size = 0.0082; % in um
xlink_cutoff = 6;

max_sites = max(n_sites);

xlinkFileStruct = '%s_protein_id.file';
%mtFileStruct = '%s_mt_coord.file';

xlinkFileName = sprintf(fileDirectory, sprintf(xlinkFileStruct, simName));
xlink_data_file = fopen(xlinkFileName);
xlink_raw_data = fread(xlink_data_file, n_mts * max_sites * n_datapoints, '*int');
fclose(xlink_data_file);
xlink_data = reshape(xlink_raw_data, max_sites, n_mts, n_datapoints);

%{
mtFileName = sprintf(fileDirectory, sprintf(mtFileStruct, simName));
mt_data_file = fopen(mtFileName);
mt_raw_data = fread(mt_data_file, n_mts * n_datapoints, '*int');
fclose(mt_data_file);
%}
mt_data = zeros(n_mts, n_datapoints); %reshape(mt_raw_data, n_mts, n_datapoints);

min_tau = time_per_datapoint;
max_tau = (n_taus - 1)*i_tau + min_tau; 
taus = min_tau:i_tau:max_tau;

dx_sq = zeros(n_taus, n_datapoints);
n_entries = zeros(n_taus, 1);

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
                    i_entry = n_entries(i_tau) + 1;
                    dx_sq(i_tau, i_entry) = dx^2;
                    n_entries(i_tau) = n_entries(i_tau) + 1;
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
                    i_entry = n_entries(i_tau) + 1;
                    dx_sq(i_tau, i_entry) = dx^2;
                    n_entries(i_tau) = n_entries(i_tau) + 1;
                else
                    disp('Greater than 2 MTs not implemented yet.');
                    return
                end
            end
        end
    end
end
%}

min_tau = time_per_datapoint;
max_tau = (n_taus - 1)*tau_increment + min_tau; 
taus = (min_tau:tau_increment:max_tau)';
MSD = zeros(n_taus, 1);
MSD_err = zeros(n_taus, 1);

for i_tau = 1:n_taus
    active_range = 1:n_entries(i_tau);
    MSD(i_tau) = mean(dx_sq(i_tau, active_range));
    variance = 0.0;
    for i_entry = active_range
        diff_sq = (MSD(i_tau) - dx_sq(i_tau, i_entry))^2;
        variance = variance + diff_sq / (n_entries(i_tau) - 1);
    end

    MSD_err(i_tau) = sqrt(variance / (n_entries(i_tau)));
end

% Use MATLAB's linear regression to find slope of data
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
