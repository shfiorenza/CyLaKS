clear all;
% Often-changed variables
simName = 'diffu_double_2_3';
n_sites = [5000, 5000];
max_sites = max(n_sites);
n_mts = length(n_sites);
starting_point = 001;
max_tau = 1.5;  % in seconds
% Pseudo-constant variables
delta_t = 0.000025;
n_steps = 40000000;
n_datapoints = 10000;
time_per_datapoint = delta_t * n_steps / n_datapoints;
active_datapoints = n_datapoints - starting_point;
site_size = 0.008; % in um
xlink_cutoff = 8;

fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
xlinkFileStruct = '%s_xlinkID.file';
mtFileStruct = '%s_mt_coord.file';

xlinkFileName = sprintf(fileDirectory, sprintf(xlinkFileStruct, simName));
mtFileName = sprintf(fileDirectory, sprintf(mtFileStruct, simName));

xlink_data_file = fopen(xlinkFileName);
xlink_raw_data = fread(xlink_data_file, n_mts * max_sites * n_datapoints, '*int');
fclose(xlink_data_file);
xlink_data = reshape(xlink_raw_data, max_sites, n_mts, n_datapoints);

mt_data_file = fopen(mtFileName);
mt_raw_data = fread(mt_data_file, n_mts * n_datapoints, '*double');
fclose(mt_data_file);
mt_data = reshape(mt_raw_data, n_mts, n_datapoints);

starting_tau = time_per_datapoint;
tau_increment = time_per_datapoint;
n_taus = max_tau/ tau_increment;

n_entries = zeros([n_taus 1]);
CSD = zeros([n_taus 1]);
MSD = zeros([n_taus 1]);

for i_tau = 1:1:n_taus
    
    tau = starting_tau + (i_tau - 1) * tau_increment;
    tau_step = int32(tau / time_per_datapoint);
    
    for i_data = (tau_step + starting_point):tau_step:n_datapoints
        % Data from current timestep
        cur_IDs = zeros(n_mts, max_sites);
        cur_coords = zeros(n_mts, 1);
        % Data for previous timestep (dictated by tau)
        prev_IDs = zeros(n_mts, max_sites);
        prev_coords = zeros(n_mts, 1);
        for i_mt=1:n_mts
            cur_IDs(i_mt, :) = xlink_data(:, i_mt, i_data);
            prev_IDs(i_mt, :) = xlink_data(:, i_mt, i_data - tau_step);
            cur_coords(i_mt) = mt_data(i_mt, i_data);
            prev_coords(i_mt) = mt_data(i_mt, i_data - tau_step);
        end
        % Scan through ID data at current timestep
        for i_site = 1:1:max_sites
            if cur_IDs(1, i_site) ~= -1
                site_coord = i_site + cur_coords(1);
                % If one MT, just search over single MT in prev. timestep
                if n_mts == 1
                    prev_index = find(prev_IDs(1, :) == cur_IDs(1, i_site));
                    if(~isempty(prev_index))
                        prev_site_coord = prev_index + prev_coords(1);
                        delta = site_coord - prev_site_coord;
                        distance = delta * site_size;
                        dist_sq = distance * distance;
                        n_entries(i_tau) = n_entries(i_tau) + 1;
                        CSD(i_tau) = CSD(i_tau) + dist_sq;
                    end
                    % Otherwise if two MTs, check to see if doubly-bound
                elseif n_mts == 2
                    i_neighb = site_coord - cur_coords(2);
                    doubly_bound = false;
                    for i_scan = -xlink_cutoff:xlink_cutoff
                        i_neighb = site_coord - cur_coords(2) + i_scan;
                        if i_neighb > 0 && i_neighb < max_sites
                            neighb_ID = cur_IDs(2, i_neighb);
                        else
                            neighb_ID = -1;
                        end
                        if neighb_ID == cur_IDs(1, i_site)
                            doubly_bound = true;
                            neighb_coord = i_neighb + cur_coords(2);
                            anchor_coord = (site_coord + neighb_coord)/2;
                            break;
                        end
                    end
                    if doubly_bound
                        prev_index = find(prev_IDs(1, :) == cur_IDs(1, i_site));
                        if(~isempty(prev_index))
                            prev_site_coord = prev_index + prev_coords(1);
                            i_neighb = prev_site_coord - prev_coords(2);
                            doubly_bound = false;
                            for i_scan = -xlink_cutoff:xlink_cutoff
                                i_neighb = prev_site_coord - prev_coords(2) + i_scan;
                                if i_neighb > 0 && i_neighb < max_sites
                                    neighb_ID = prev_IDs(2, i_neighb);
                                else
                                    neighb_ID = -1;
                                end
                                if neighb_ID == prev_IDs(1, i_site)
                                    doubly_bound = true;
                                    prev_neighb_coord = i_neighb + prev_coords(2);
                                    prev_anchor_coord = (prev_site_coord + prev_neighb_coord)/2;
                                    break;
                                end
                            end
                            if doubly_bound
                                delta = anchor_coord - prev_anchor_coord;
                                %msg = sprintf('delta is %f for xlink %i | tau_step %i | i_data %i', ...
                                %    delta, cur_IDs(1, i_site), tau_step, i_data);
                                %disp(msg);
                                distance = delta * site_size;
                                dist_sq = distance * distance;
                                n_entries(i_tau) = n_entries(i_tau) + 1;
                                CSD(i_tau) = CSD(i_tau) + dist_sq;
                            end
                            
                        end
                    end
                else
                    disp('Greater than 2 MTs not implemented yet.');
                end
            end
        end
    end
    %loc_entries = n_entries(i_tau);
    %loc_CSD = CSD(i_tau);
    MSD(i_tau) = CSD(i_tau) / n_entries(i_tau);
end

% Use basic linear regression to find slope
y = MSD;
x = zeros(length(MSD), 1);
for i_tau = 1:1:n_taus
    x(i_tau) = starting_tau + (i_tau - 1) * tau_increment;
end

X = [ones(length(x), 1) x];
m = X\y;
D = m(2) / 2

%d_consts(i_eff) = D;

fig1 = figure();
set(fig1, 'Position', [50, 50, 960, 600])

y2 = X*m;
plot(x, y2, '--', 'LineWidth',  2)
hold on
plot(x, y, '*', 'LineWidth', 2)
title({'Two Anti-Parallel Microtubules','c_{bulk} = 0.05 nM  |  L_{MTs} = 400 \mum  |  D_{ii} = 0.0437 \mum^2s^{-1}'});
ylabel('MSD (\mum^2)');
xlabel('Tau (s)');
xlim([0.0 (max_tau + starting_tau)]);
legend('Fit', 'Data', 'location', 'northeastoutside');

dim = [0.18 0.75 .1 .1];
str = sprintf('D_{obs} = %#.3g um^2s^{-1}', D);
annotation('textbox',dim,'String',str,'FitBoxToText','on');

%end
%{
fig1 = figure();
set(fig1, 'Position', [50, 50, 960, 600]);
semilogx(conc_eff, d_consts);
hold on
line([0, conc_eff(length(conc_eff))], [0.0109, 0.0109]);
title('c = 0.05 nM;  D\_i = 0.038 um^2/;  D\_ii = 0.019 um^2/s');
ylabel('Diffusion constant (um^2)/s');
xlabel('Effective binding concentration (nM)');
grid on
grid minor
%}