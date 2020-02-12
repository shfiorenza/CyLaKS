
clear all;
simNameBase = "mt_diffusion";
mt_lengths = [125, 500, 2000]; % in n_sites
D_expected = [1.36, 0.340, 0.0851];
%{
mt_lengths = [125, 250, 500, 1000, 2000, 4000]; % in n_sites
D_expected = [1.36, 0.681, 0.340, 0.170, 0.0851, 0.0423];
%}
seeds = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
n_mts = length(mt_lengths);
n_seeds = length(seeds);
max_tau = 48.015;
% Pseudo-constant variables
site_size = 0.008; % in um
delta_t = 2.5e-5;
n_steps = 60000000;
n_datapoints = 100000;
time_per_datapoint = delta_t * n_steps / n_datapoints;

fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
mtFileStruct = '%s_mt_coord.file';

min_tau = time_per_datapoint;
tau_increment = 200*time_per_datapoint;
n_taus = (max_tau - min_tau) / tau_increment;
max_tau_step = max_tau / time_per_datapoint;


squared_displacements = zeros(n_mts, n_taus, n_datapoints * n_seeds);
n_entries = zeros(n_mts, n_taus);
for i_seed = 1 : n_seeds
    simName = simNameBase;
    if length(seeds) > 1
        simName = sprintf("%s_%i", simNameBase, i_seed - 1);
    end
    mtFileName = sprintf(fileDirectory, sprintf(mtFileStruct, simName));
    mt_data_file = fopen(mtFileName);
    mt_raw_data = fread(mt_data_file, n_mts * n_datapoints, '*int');
    fclose(mt_data_file);
    mt_data = reshape(mt_raw_data, n_mts, n_datapoints);
    for i_mt = 1 : n_mts
        for i_tau = 1 : n_taus
            tau = min_tau + (i_tau - 1) * tau_increment;
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
MSD = zeros(n_mts, n_taus);
MSD_err = zeros(n_mts, n_taus);
for i_mt = 1 : n_mts
    for i_tau = 1 : n_taus
        active_range = 1:n_entries(i_mt, i_tau);
        MSD(i_mt, i_tau) = mean(squared_displacements(i_mt, i_tau, active_range));
        variance = 0.0;
        for i_entry = active_range
            diff_sq = (MSD(i_mt, i_tau) - squared_displacements(i_mt, i_tau, i_entry))^2;
            variance = variance + diff_sq;
        end
        variance = variance/double(n_entries(i_mt, i_tau) - 1);
        MSD_err(i_mt, i_tau) = sqrt(variance);
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

fig1 = figure();
set(fig1, 'Position', [50, 50, 1500, 400]);
taus = linspace(min_tau, max_tau, n_taus);

for i_mt = 1 : n_mts
subplot(1, n_mts, i_mt)
errorbar(taus, MSD(i_mt, :), MSD_err(i_mt, :), 'o', 'LineWidth', 2); %, 'Color', [0 0.447 0.741])
hold on
plot([0 max_tau], [0 2*D_expected(i_mt)*max_tau], '--', 'LineWidth', 2); %, 'Color', [0.85 0.325 0.098]); 
ylabel("Mean squared displacement (um^2)");
xlabel("Tau (s)");
if mt_lengths(i_mt) * site_size == 1
    title(sprintf("Microtubule length = %g micron", mt_lengths(i_mt) * site_size));
else
    title(sprintf("Microtubule length = %g microns", mt_lengths(i_mt) * site_size));
end
legend(["Simulation data", "Theory"], 'location', 'northwest');
end
%{
subplot(1, 3, 2)
errorbar(taus, MSD(3, :), MSD_err(3, :), 'o', 'LineWidth', 2, 'Color', [0.4940    0.1840    0.5560])
hold on
plot([0 max_tau], [0 2*D_expected(3)*max_tau], '--', 'LineWidth', 2, 'Color', [0.4660    0.6740    0.1880]); 
ylabel("Mean squared displacement (um^2)");
xlabel("Tau (s)");
title(sprintf("Length = %g microns", mt_lengths(3) * site_size));
legend(["Simulation data", "Theory"], 'location', 'northwest');

subplot(1, 3, 3)
errorbar(taus, MSD(6, :), MSD_err(6, :), 'o', 'LineWidth', 2, 'Color', [0.6350    0.0780    0.1840])
hold on
plot([0 max_tau], [0 2*D_expected(6)*max_tau], '--', 'LineWidth', 2, 'Color', [0.3010    0.7450    0.9330]); 
ylabel("Mean squared displacement (um^2)");
xlabel("Tau (s)");
title(sprintf("Length = %g microns", mt_lengths(6) * site_size));
legend(["Simulation data", "Theory"], 'location', 'northwest');
%}
%{
for i_mt = 1 : 1 : n_mts
    %errorbar(taus, MSD(i_mt, :), MSD_err(i_mt, :), 'LineWidth', 2, 'Color', RGBArray(i_mt,:));
    plot(taus, MSD(i_mt, :), 'LineWidth', 2, 'Color', RGBArray(i_mt, :));
end
legendLabel = "Length = " + num2str(mt_lengths'* site_size, "%i") + " microns";
legend(legendLabel, 'location', 'northwest', 'AutoUpdate', 'off');
for i_mt = 1 : 1 : n_mts
    plot([0 max_tau], [0 2*D_expected(i_mt)*max_tau], '--', 'LineWidth', 2, ...
        'Color', RGBArray(i_mt, :)); 
end
%}
