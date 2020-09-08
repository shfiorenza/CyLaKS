
baseNames = ["lattice_coop"];
concentrations = [20, 50, 80, 120, 220, 420];
dir = "/home/shane/Projects/overlap_analysis/mgh_model/";

n_runs = length(baseNames);
n_concs = length(concentrations);

runlengths = zeros(n_runs, n_concs);
err_runlengths = zeros(n_runs, n_concs);
lifetimes = zeros(n_runs, n_concs);
err_lifetimes = zeros(n_runs, n_concs);
velocities = zeros(n_runs, n_concs);
err_velocities = zeros(n_runs, n_concs);

for i_run = 1 : n_runs
    for i_conc = 1 : n_concs
        conc = int32(concentrations(i_conc));
        simName = sprintf("%s/%s_%i", dir, baseNames(i_run), conc);
        mot_stats = get_motor_stats(simName);
        runlengths(i_run, i_conc) = mot_stats(1);
        err_runlengths(i_run, i_conc) = mot_stats(2);
        lifetimes(i_run, i_conc) = mot_stats(3);
        err_lifetimes(i_run, i_conc) = mot_stats(4);
        velocities(i_run, i_conc) = mot_stats(5);
        err_velocities(i_run, i_conc) = mot_stats(6);
    end
end

%{
% best_wtSq_ATP; stats: [0.76, 851, 0.75, 3.75] (0.175, 2.3) 
runlengths = [1.37, 1.40, 1.63, 2.62, 1.97, 1.30];
err_runlengths = [0.2, 0.16, 0.14, 0.16, 0.1, 0.08];
lifetimes = [2.33, 2.37, 3.23, 13.6, 11.4, 10.7];
err_lifetimes = [0.4, 0.3, 0.3, 0.8, 0.7, 0.7];
velocities = [590, 590, 580, 440, 320, 190];
err_velocities = [12, 8, 10, 15, 15, 13];

%best_allSix; stats: [1.12, 872, 0.504, 4.08, 0.124, 3.16]
runlengthsB = [1.15, 1.17, 1.36, 1.46, 1.99, 1.21];
err_runlengthsB = [0.2, 0.1, 0.1, 0.1, 0.1, 0.08];
lifetimesB = [2.02, 2.05, 2.42, 2.58, 11.1, 10.5];
err_lifetimesB = [0.4, 0.2, 0.2, 0.2, 0.7, 0.7];
velocitiesB = [570, 570, 570, 560, 390, 230];
err_velocitiesB = [16, 10, 10, 10, 15, 15];
%}

%{
% best_wtSq; stats: [2.8, 434, 1.03, 3.4]
runlengths = [0.956, 1.04, 1.05, 1.11, 2.66, 3.07];
err_runlengths = [0.17, 0.13, 0.12, 0.12, 0.26, 0.25];
lifetimes = [1.44, 1.55, 1.65, 1.65, 12.3, 26.5];
err_lifetimes = [0.25, 0.2, 0.19, 0.18, 1.2, 2.1];
velocities = [665, 663, 676, 673, 517, 350];
err_velocities = [18, 10, 12, 10, 25, 24];

%best; stats: [2.8, 875, 0.658, 3.6]
runlengthsB = [0.97, 1.03, 0.96, 1.03, 1.19, 7.55];
err_runlengthsB = [0.17, 0.13, 0.11, 0.11, 0.13, 0.56];
lifetimesB = [1.45, 1.54, 1.45, 1.53, 1.74, 19.05];
err_lifetimesB = [0.25, 0.2, 0.17, 0.17, 0.19, 1.4];
velocitiesB = [674, 674, 670, 681, 683, 526];
err_velocitiesB = [15, 11, 12, 11, 10, 14];
%}
%{
% velocity-only fit
runlengths = [0.9,  1.1, 3.4, 4.6, 5.0, 4.5];
err_runlengths = [0.17, 0.15, 0.4, 0.5, 0.6, 0.5];
lifetimes = [1.5, 1.8, 10.4, 18.3, 26.3, 35.5];
err_lifetimes = [0.3, 0.2, 1, 2, 3, 4];
velocities = [621, 616, 560, 451, 327, 196];
err_velocities = [10, 10, 20, 20, 20, 10];

%lifetime-only fit
runlengthsB = [0.9, 0.9, 1.0, 1.0, 2.2, 3.9];
err_runlengthsB = [0.2, 0.1, 0.14, 0.1, 0.3, 0.5];
lifetimesB = [1.5, 1.5, 1.7, 1.6, 4.7, 20.9];
err_lifetimesB = [0.3, 0.2, 0.2, 0.2, 0.7, 2.7];
velocitiesB = [622, 624, 620, 620, 561, 396];
err_velocitiesB = [10, 10, 10, 10, 20, 30];
%}
%{
% without cap
runlengthsB = [0.9, 1.0 ,1.0, 1.5, 7.1, 5.8];
err_runlengthsB = [0.08, 0.07, 0.06, 0.09, 0.3, 0.3];
lifetimesB = [1.5, 1.6, 1.6, 2.7, 28.2, 40.5];
err_lifetimesB = [0.1, 0.1, 0.09, 0.2, 1, 2];
velocitiesB = [628, 625, 624, 610, 370, 280];
err_velocitiesB = [6, 4, 4, 5, 7, 7];

% with cap
runlengths = [0.9, 1.0 ,1.0, 1.5, 6.8, 5.1];
err_runlengths = [0.08, 0.07, 0.06, 0.09, 0.3, 0.3];
lifetimes = [1.5, 1.6, 1.6, 2.7, 18.8, 16.2];
err_lifetimes = [0.1, 0.1, 0.09, 0.2, 0.9, 1];
velocities = [628, 625, 624, 610, 410, 350];
err_velocities = [6, 4, 4, 5, 6, 6];
%}

exp_concs = [20, 50, 80, 120, 220, 420];
exp_runlengths = [0.970, 1.310, 2.420, 1.660, 1.960, 2.86];
exp_err_runlengths = [0.180, 0.320, 0.350, 0.940, 0.310, 0.72];
exp_lifetimes = [1.8, 2.1, 7.1, 5.2, 8.3, 17.9];
exp_err_lifetimes = [0.6, 0.7, 1.7, 5.9, 2.6, 3.9];
exp_velocities = [600, 710, 360, 310, 310, 180];
exp_err_velocities = [75, 110, 50, 79, 40, 40];

fig1 = figure();
set(fig1, 'Position', [50, 50, 2*720, 2*240])
    
subplot(1, 4, 1)
errorbar(exp_concs, exp_runlengths, exp_err_runlengths, ':d','LineWidth', 2);
hold on
for i_run = 1 : n_runs
    errorbar(concentrations, runlengths(i_run, :), err_runlengths(i_run, :), 'o','LineWidth', 2);
end
%errorbar(concentrations, runlengthsB, err_runlengthsB, 'o','LineWidth', 2);
ylabel('Run length (microns)');

subplot(1, 4, 2)
errorbar(exp_concs, exp_lifetimes, exp_err_lifetimes, ':d','LineWidth', 2);
hold on
for i_run = 1 : n_runs
    errorbar(concentrations, lifetimes(i_run, :), err_lifetimes(i_run, :), 'o', 'LineWidth', 2);
end
%errorbar(concentrations, lifetimesB, err_lifetimesB, 'o', 'LineWidth', 2);
xlabel('KIF4A concentration (pM)');
ylabel('Life time (seconds)');

subplot(1, 4, 3)
errorbar(exp_concs, exp_velocities, exp_err_velocities, ':d','LineWidth', 2);
hold on
for i_run = 1 : n_runs
    errorbar(concentrations, velocities(i_run, :), err_velocities(i_run, :), 'o', 'LineWidth', 2);
end
%errorbar(concentrations, velocitiesB, err_velocitiesB, 'o', 'LineWidth', 2);
ylabel('Velocity (nm/s)');

subplot(1, 4, 4);
hold on
errorbar([0, 0], [0, 0], '--o', 'LineWidth', 2);
errorbar([0, 0], [0, 0], '--o', 'LineWidth', 2);
%errorbar([0, 0], [0, 0], '--o', 'LineWidth', 2);
xlim([1 2]);
ylim([1 2]);
axis off
legend({'Data', 'Best fit (4 params)', 'Best fit (6 params)'}, 'location', 'northwestoutside'); %,  ...
    %'Sim data w/ cooperative unbind\_ii'}, 'location', 'best');
