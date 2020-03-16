
clear all;
% Dynamic variables
simNameBase = "mt_diffusion";
%mt_lengths = [25, 125, 500, 2000, 4000]; % in n_sites
%D_expected = [6.80771, 2.41534, 0.83076, 0.264421, 0.146393];
mt_lengths = [25, 125, 625, 3125, 15625];
D_expected = [6.80771, 2.41534, 0.693829, 0.180918, 0.044614];
seeds = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
tau_increment = 9;
max_tau = 99.15;
% Pseudo-constant variables
site_size = 0.008;      % microns
delta_t = 2.5e-5;       % seconds
n_steps = 60000000;
n_datapoints = 10000;
% Quantities calculated from input variables
time_per_datapoint = delta_t * n_steps / n_datapoints;
min_tau = time_per_datapoint;
taus = min_tau : tau_increment : max_tau;
n_taus = length(taus);
n_seeds = length(seeds);
n_mts = length(mt_lengths);
% File directory stuff
fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
mtFileStruct = '%s_mt_coord.file';

% Run thru mt_coord data and get squared displacements for all taus
squared_displacements = zeros(n_mts, n_taus, n_datapoints * n_seeds);
n_entries = zeros(n_mts, n_taus);
for i_seed = 1 : n_seeds
    simName = simNameBase;
    if length(seeds) > 1
        simName = sprintf("%s_%i", simNameBase, seeds(i_seed));
    end
    mtFileName = sprintf(fileDirectory, sprintf(mtFileStruct, simName));
    mt_data_file = fopen(mtFileName);
    mt_data = fread(mt_data_file, n_mts * n_datapoints, '*int');
    fclose(mt_data_file);
    mt_data = reshape(mt_data, n_mts, n_datapoints);
    for i_mt = 1 : n_mts
        for i_tau = 1 : n_taus
            tau = taus(i_tau);
            tau_step = int32(tau / time_per_datapoint);  
            for i_data = 1 : (n_datapoints - tau_step)
                cur_pos = mt_data(i_mt, i_data);
                next_pos = mt_data(i_mt, i_data + tau_step);
                dx = double(next_pos - cur_pos) * site_size;
                i_entry = n_entries(i_mt, i_tau) + 1;
                squared_displacements(i_mt, i_tau, i_entry) = dx^2;
                n_entries(i_mt, i_tau) = n_entries(i_mt, i_tau) + 1;
            end
        end
    end
end
% Calculate MSD and standard error of MSD from squared displacement entries
MSD = zeros(n_mts, n_taus);
MSD_err = zeros(n_mts, n_taus);
for i_mt = 1 : n_mts
    for i_tau = 1 : n_taus
        active_range = 1:n_entries(i_mt, i_tau);
        MSD(i_mt, i_tau) = mean(squared_displacements(i_mt, i_tau, active_range));
        variance = 0.0;
        for i_entry = active_range
            diff_sq = (MSD(i_mt, i_tau) - squared_displacements(i_mt, i_tau, i_entry))^2;
            variance = variance + diff_sq / double(n_entries(i_mt, i_tau) - 1);
        end
        MSD_err(i_mt, i_tau) = sqrt(variance / double(n_entries(i_mt, i_tau)));
    end
end

%{
D = zeros(n_mts, 1);
y_int = zeros(n_mts, 1);
for i_mt=1:n_mts
    % Use basic linear regression to find slope
    y = MSD(i_mt, :)';
    x = zeros(length(MSD(i_mt, :)), 1);
    for i_tau = 1:1:n_taus
        x(i_tau) = min_tau + (i_tau - 1) * tau_increment;
    end
    X = [ones(length(x), 1) x];
    m = X\y;
    D(i_mt) = m(2) / 2;
    y_int(i_mt) = m(1);
    fprintf("MT %i: D = %g, y_int = %g\n", i_mt, D(i_mt), y_int(i_mt));
end
%}
%}

% Plot
fig1 = figure();
set(fig1, 'Position', [50, 50, 1500, 500]);
for i_mt = 1 : n_mts
    subplot(1, n_mts + 1, i_mt)
    hold on
    errorbar(taus, MSD(i_mt, :), MSD_err(i_mt, :), 'o', 'LineWidth', 2, 'MarkerSize', 10);
    plot([0 max_tau], [0 2*D_expected(i_mt)*max_tau], '--', 'LineWidth', 2); %, ...
       % 'Color', [0.5 0.5 0.5]);
    ax = gca;
    ax.FontSize = 10;   
    xlabel("Tau (s)", 'FontSize', 12);
    ylabel("Mean squared displacement (\mum^2)", 'FontSize', 12);
    title({sprintf("Length = %g", mt_lengths(i_mt) * site_size) + " \mum", " "}, 'FontSize', 12);
    
end
% Some trickery to give all plots a common legend
subplot(1, n_mts + 1, n_mts + 1)
hold on
errorbar([0,0], [0, 0], 'o', 'LineWidth', 2, 'MarkerSize', 10);
plot([0, 0], '--', 'LineWidth', 2);
xlim([1 2]);
ylim([1 2]);
axis off
legend(["Sim data", "Theory"], 'location', 'northwestoutside', 'FontSize', 10);