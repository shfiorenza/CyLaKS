%{
clear variables;
file_dir = '/home/shane/projects/CyLaKS/%s';
sim_name_base = 'run_mt_diffusion/mt_diffusion';
seeds = [0, 1, 2, 3, 4, 5];
n_taus = 10;
i_tau = 10; % 10.0;

n_dims = 2;
% Open log file and parse it into param labels & their values
log_file = sprintf(file_dir, sprintf('%s.log', sim_name_base));
if(~isempty(seeds))
    log_file = sprintf(file_dir, sprintf('%s_%i.log', sim_name_base, seeds(1)));
end
log = textscan(fileread(log_file), '%s %s', 'Delimiter', '=');
params = log{1, 1};
values = log{1, 2};
% Read in number of MTs
n_mts = sscanf(values{contains(params, 'count ')}, '%g');
n_sites = zeros(1, n_mts);
ell = zeros(1, n_mts);
D_par = zeros(1, n_mts);
D_perp = zeros(1, n_mts);
for i_mt = 1 : n_mts
    string = sprintf('n_sites[%i] ', i_mt - 1);
    n_sites(i_mt) = sscanf(values{contains(params, string)}, '%i');
    string = sprintf('length[%i] ', i_mt - 1);
    ell(i_mt) = sscanf(values{contains(params, string)}, '%g');
    string = sprintf('D_par[%i] ', i_mt - 1);
    D_par(i_mt) = sscanf(values{contains(params, string)}, '%g');
    string = sprintf('D_perp[%i] ', i_mt - 1);
    D_perp(i_mt) = sscanf(values{contains(params, string)}, '%g');
end
% Read in system params
dt = sscanf(values{contains(params, 'dt ')}, '%g');
steps_per_datapoint = str2double(values{contains(params, 'n_steps_per_snapshot ')});
time_per_datapoint = dt * steps_per_datapoint;
n_datapoints = str2double(values{contains(params, 'n_datapoints ')});
% Use actual recorded number of datapoints to parse thru data/etc
if any(contains(params, 'N_DATAPOINTS ') ~= 0)
    n_datapoints = str2double(values{contains(params, 'N_DATAPOINTS ')});
end
n_seeds = 1;
if(~isempty(seeds))
    n_seeds = length(seeds);
end
min_tau = time_per_datapoint;
max_tau = (n_taus - 1)*i_tau + min_tau; 
taus = min_tau:i_tau:max_tau;

com_x = zeros(n_mts, n_datapoints);
com_y = zeros(n_mts, n_datapoints);
% Run thru mt_coord data and get squared displacements for all taus
dx_sq = zeros(n_mts, n_taus, n_datapoints * n_seeds);
dy_sq = zeros(n_mts, n_taus, n_datapoints * n_seeds);
n_entries = zeros(n_mts, n_taus);
% Scan thru different seeds
for i_seed = 1 : n_seeds
    sim_name = sim_name_base; 
    if(~isempty(seeds))
        sim_name = sprintf('%s_%i', sim_name_base, seeds(i_seed));
    end
    filename = '%s_filament_pos.file';
    file = fopen(sprintf(file_dir, sprintf(filename, sim_name)));
    data = fread(file, 2*n_dims * n_mts * n_datapoints, '*double');
    filament_pos = reshape(data, n_dims, 2, n_mts, n_datapoints);
    fclose(file);
    % Scan over individual MTs
    for i_mt = 1 : n_mts
        for i_data = 1 : n_datapoints
            plus_pos = filament_pos(:, 1, i_mt, i_data);
            minus_pos = filament_pos(:, 2, i_mt, i_data);
            com_x(i_mt, i_data) = double(plus_pos(1) + minus_pos(1))/2.0;
            com_y(i_mt, i_data) = double(plus_pos(2) + minus_pos(2))/2.0;
        end
        % Scan over tau values
        for i_tau = 1 : n_taus
            tau = taus(i_tau);
            tau_step = int32(tau / time_per_datapoint);
            for i_data = 1 : (n_datapoints - tau_step)
                plus_pos = filament_pos(:, 1, i_mt, i_data);
                minus_pos = filament_pos(:, 2, i_mt, i_data);
                x = com_x(i_mt, i_data);
                y = com_y(i_mt, i_data);
                plus_pos_tau = filament_pos(:, 1, i_mt, i_data + tau_step);
                minus_pos_tau = filament_pos(:, 2, i_mt, i_data + tau_step);
                x_tau = com_x(i_mt, i_data + tau_step); %(plus_pos_tau(1) + minus_pos_tau(1))/2;
                y_tau = com_y(i_mt, i_data + tau_step);%(plus_pos_tau(2) + minus_pos_tau(2))/2;
                
                i_entry = n_entries(i_mt, i_tau) + 1;
                dx_sq(i_mt, i_tau, i_entry) = (x - x_tau)^2;
                dy_sq(i_mt, i_tau, i_entry) = (y - y_tau)^2;
                n_entries(i_mt, i_tau) = n_entries(i_mt, i_tau) + 1;
            end
        end
    end
end

% Calculate MSD and standard error of MSD from squared displacement entries
MSD_par = zeros(n_mts, n_taus);
MSD_perp = zeros(n_mts, n_taus);
MSD_par_err = zeros(n_mts, n_taus);
MSD_perp_err = zeros(n_mts, n_taus);
for i_mt = 1:n_mts
    for i_tau = 1:n_taus
        active_range = 1:n_entries(i_mt, i_tau);
        MSD_par(i_mt, i_tau) = mean(dx_sq(i_mt, i_tau, active_range));
        MSD_perp(i_mt, i_tau) = mean(dy_sq(i_mt, i_tau, active_range));
        var_par = 0.0;
        var_perp = 0.0;
        for i_entry = active_range
            diff_sq = (MSD_par(i_mt, i_tau) - dx_sq(i_mt, i_tau, i_entry))^2;
            var_par = var_par + diff_sq / double(n_entries(i_mt, i_tau) - 1);
            diff_sq = (MSD_perp(i_mt, i_tau) - dy_sq(i_mt, i_tau, i_entry))^2;
            var_perp = var_perp + diff_sq / double(n_entries(i_mt, i_tau) - 1);
        end
        MSD_par_err(i_mt, i_tau) = sqrt(var_par / double(n_entries(i_mt, i_tau)));
        MSD_perp_err(i_mt, i_tau) = sqrt(var_perp / double(n_entries(i_mt, i_tau)));
    end
end

%{
D = zeros(n_mts, 1);
y_int = zeros(n_mts, 1);

for i_mt = 1:n_mts
    % Use basic linear regression to find slope
    y = MSD(i_mt, :)';
    x = zeros(length(MSD(i_mt, :)), 1);

    for i_tau = 1:1:n_taus
        x(i_tau) = min_tau + (i_tau - 1) * tau_increment;
    end

    X = [ones(length(x), 1) x];
    m = X \ y;
    D(i_mt) = m(2) / 2;
    y_int(i_mt) = m(1);
    fprintf("MT %i: D = %g, y_int = %g\n", i_mt, D(i_mt), y_int(i_mt));
end
%}
%}

% Plot
fig1 = figure();
set(fig1, 'Position', [50, 50, 720, 540]);

data_color = [0, 0.4470, 0.7410; 0.9290, 0.6940, 0.1250; 0.4660, 0.6740, 0.1880];
line_color = [0.8500, 0.3250, 0.0980; 0.4940, 0.1840, 0.5560; 0.3010, 0.7450, 0.9330];
color = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980;  0.9290, 0.6940, 0.1250];

c = 1e-6; % convert from nm^2 to um^2

for i_mt = 1:n_mts
   % subplot(1, n_mts, i_mt)
    hold on
    par_data = plot(taus, c*MSD_par(i_mt, :), 'o', 'MarkerSize', 12, 'MarkerEdgeColor', color(i_mt, :));
    par_data.MarkerFaceColor = par_data.MarkerEdgeColor;
    par_data.Color = par_data.MarkerFaceColor;
    %errorbar(taus, MSD_par(i_mt, :), MSD_par_err(i_mt, :), 'o', 'LineWidth', 2, 'MarkerSize', 10);
    perp_data = plot(taus, c*MSD_perp(i_mt, :), 'sq', 'MarkerSize', 12, 'MarkerEdgeColor', color(i_mt, :));
    perp_data.MarkerFaceColor = perp_data.MarkerEdgeColor;
    perp_data.Color = perp_data.MarkerFaceColor;
   plot([0 max_tau + i_tau/8], [0 c*2 * D_par(i_mt) * (max_tau + i_tau/8)], '-', 'LineWidth', 2, 'Color', [0.6 0.6 0.6]); %line_color(i_mt, :)); 
     plot([0 max_tau + i_tau/8], [0 c*2 * D_perp(i_mt) * (max_tau + i_tau/8)], '-', 'LineWidth', 2, 'Color', [0.6 0.6 0.6]); %line_color(i_mt, :));
    
    %errorbar(taus, MSD_perp(i_mt, :), MSD_perp_err(i_mt, :), 'o', 'LineWidth', 2, 'MarkerSize', 10);
    set(gca, 'FontSize', 22);
    xlabel("Tau (s)", 'FontSize', 22);
    ylabel("MSD (\mum^2)", 'FontSize', 22);
    xlim([-i_tau/4 max_tau + i_tau/4]);
    ylim([-0.15 5]);
    xticks([0 25 50 75]);
    yticks([0 2 4]);
end
h = get(gca,'Children');
n_entries = length(h);
h_array = [];
for i = n_entries : - 1 : 1
      h_array = [h_array h(i)];
end
set(gca,'Children',h_array)
legend(h([n_entries, n_entries-4, n_entries-8, n_entries-2, n_entries-1, n_entries-5, n_entries-9]), ... 
    ["L = 1 um (par)", "L = 5 um (par)", "L = 20 um (par)", "Theory", "L = 1um (perp)", "L = 5 um (perp)", "L = 20 um (perp)"], ...
    'location', 'northwest', 'FontSize', 18, 'NumColumns', 2);
legend('boxoff') 
%{
% Some trickery to give all plots a common legend
subplot(1, n_mts + 1, n_mts + 1)
hold on
errorbar([0, 0], [0, 0], 'o', 'LineWidth', 2, 'MarkerSize', 10);
plot([0, 0], '--', 'LineWidth', 2);
errorbar([0, 0], [0, 0], 'o', 'LineWidth', 2, 'MarkerSize', 10);
plot([0, 0], '--', 'LineWidth', 2);
xlim([1 2]);
ylim([1 2]);
axis off

%}
