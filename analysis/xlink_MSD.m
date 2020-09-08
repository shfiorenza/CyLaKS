clear all;
% Often-changed variables
n_sites = 500;
%simName = 'test6c';
% Pseudo-constant variables
delta_t = 0.0005;
n_steps = 2000000;
n_mts = 2;
n_datapoints = 100000;
time_per_datapoint = delta_t * n_steps / n_datapoints;
starting_point = 1;
active_datapoints = n_datapoints - starting_point;
site_size = 0.008; % in um
xlink_cutoff = 7; 

fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
xlinkFileStruct = '%s_xlinkID.file';
mtFileStruct = '%s_mt_coord.file';

conc_eff = [25,50,250,500,2500,5000,25000,50000,250000,500000,2500000,5000000];

d_consts = zeros([length(conc_eff) 1]);

for i_eff=1:1:length(conc_eff)

simName = sprintf('xl_scan/XlinkDiffScan_%i',conc_eff(i_eff));
    
xlinkFileName = sprintf(fileDirectory, sprintf(xlinkFileStruct, simName));
mtFileName = sprintf(fileDirectory, sprintf(mtFileStruct, simName));

xlink_data_file = fopen(xlinkFileName);
xlink_raw_data = fread(xlink_data_file, [n_mts * n_sites * n_datapoints], '*int');
fclose(xlink_data_file);
xlink_data = reshape(xlink_raw_data, n_sites, n_mts, n_datapoints);

%mt_data_file = fopen(mtFileName);
%mt_raw_data = fread(mt_data_file, [n_mts * n_datapoints], '*double');
%fclose(mt_data_file);
mt_data = zeros([n_mts n_datapoints]); %reshape(mt_raw_data, n_mts, n_datapoints);

max_tau = 4;  % in seconds
n_taus = 40;
tau_increment = max_tau / n_taus;
starting_tau = 3;

n_entries = zeros([n_taus 1]);
CSD = zeros([n_taus 1]); 
MSD = zeros([n_taus 1]);

for i_tau = starting_tau:1:n_taus
    
    tau = i_tau * tau_increment;
    tau_step = int32(tau / time_per_datapoint);
    
    for i_data = starting_point:tau_step:active_datapoints
        % Ensure enough time has passed to have 2 data points
        if(i_data > tau_step)
            % Import data for current timestep
            cur_IDs_mtOne = xlink_data(:, 1, i_data);
            cur_IDs_mtTwo = xlink_data(:, 2, i_data);
            cur_mtOne_coord = mt_data(1, i_data);
            cur_mtTwo_coord = mt_data(2, i_data); 
            % Import data for previous timestep (dictated by tau)
            prev_IDs_mtOne = xlink_data(:, 1, i_data - tau_step);
            prev_IDs_mtTwo = xlink_data(:, 2, i_data - tau_step); 
            prev_mtOne_coord = mt_data(1, i_data - tau_step);
            prev_mtTwo_coord = mt_data(2, i_data - tau_step); 
            % Scan through ID data at current timestep
            for i_site = 1:1:n_sites 
                xlink_ID = cur_IDs_mtOne(i_site);
                % Make sure an xlink occupies this site
                if(xlink_ID ~= -1)
                    % Make sure xlink is doubly-bound
                    site_coord = i_site + cur_mtOne_coord;
                    neighb_index = site_coord - cur_mtTwo_coord;
                    double_bound = false; 
                    cur_anchor_pos = 0; 
                    for i_scan = -xlink_cutoff:1:xlink_cutoff
                        scan_index = neighb_index + i_scan;
                        if(scan_index > 0 && scan_index < n_sites)
                            neighb_ID = cur_IDs_mtTwo(neighb_index + i_scan);
                        else
                            neighb_ID = -1;
                        end
                        % If doubly-bound, find current anchor coordinate
                        if(neighb_ID == xlink_ID)
                            double_bound = true;
                            neighb_coord = neighb_index + i_scan + cur_mtTwo_coord;
                            cur_anchor_pos = (site_coord + neighb_coord)/2;
                            break;
                        end
                    end 
                    % If xlink was doubly-bound, search for xlink ID at previous timestep 
                    if(double_bound == true)
                        prev_index = find(prev_IDs_mtOne == xlink_ID);
                        % If found, repeat above process
                        if(~isempty(prev_index))
                            site_coord = prev_index + prev_mtOne_coord;
                            neighb_index = site_coord - prev_mtTwo_coord;
                            double_bound = false;
                            prev_anchor_pos = 0;
                            for i_scan = -xlink_cutoff:1:xlink_cutoff
                                scan_index = neighb_index + i_scan;
                                if(scan_index > 0 && scan_index < n_sites)
                                    neighb_ID = prev_IDs_mtTwo(neighb_index + i_scan);
                                else
                                    neighb_ID = -1;
                                end
                                % If doubly-bound, find prev anchor coord
                                if(neighb_ID == xlink_ID)
                                    double_bound = true;
                                    neighb_coord = neighb_index + i_scan + prev_mtTwo_coord;
                                    prev_anchor_pos = (site_coord + neighb_coord)/2;
                                    break;
                                end
                            end
                            if(double_bound == true)
                                delta = cur_anchor_pos - prev_anchor_pos;
                                distance = delta * site_size;
                                dist_sq = distance * distance;
                                n_entries(i_tau) = n_entries(i_tau) + 1;
                                CSD(i_tau) = CSD(i_tau) + dist_sq;
                            end
                        end
                    end
                end
            end
        end
    end
    loc_entries = n_entries(i_tau);
    loc_CSD = CSD(i_tau);
    MSD(i_tau) = CSD(i_tau) / n_entries(i_tau);
end

% Use basic linear regression to find slope
y = MSD(starting_tau:n_taus);
x = zeros([length(y) 1]);

for i_tau = starting_tau:1:n_taus
    x(i_tau-starting_tau + 1) = i_tau * tau_increment;
end

X = [ones(length(x), 1) x];
m = X\y;

D = m(2) / 2;

d_consts(i_eff) = D;
%{
fig1 = figure();
set(fig1, 'Position', [50, 50, 960, 600])

y2 = X*m;
plot(x, y2, 'LineWidth', 2)
hold on
plot(x, y, 'LineWidth', 2)
title('c = 0.05 nM;  eff\_bind\_conc = 750 nM;  D\_i = 0.036 um^2/;  D\_ii = 0.018 um^2/s');
ylabel('MSD (um^2)');
xlabel('Tau (s)');
legend('Fit', 'Data', 'location', 'northeastoutside');

dim = [0.21 0.75 .1 .1];
str = sprintf('D_{eff} = %#.3g um^2/s', D);
annotation('textbox',dim,'String',str,'FitBoxToText','on');
%}
end

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