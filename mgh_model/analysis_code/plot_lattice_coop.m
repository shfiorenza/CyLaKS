
clear all;

baseNames = ["lattice_coop_base"]; %, "lattice_coop_justBind"]; %, "lattice_coopH"];
%baseNames = ["lattice_coopH", "lattice_coopM", "lattice_coopN"];
labelNames = ["Sim. w/ short-range only", "Sim. w/ long-range", ...
              "Sim. w/ long-range (also affects stepping)"];
    
concentrations = [20, 50, 80, 120, 220, 420];
dir = "/home/shane/Projects/overlap_analysis/mgh_model";

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

exp_concs = [20, 50, 80, 120, 220, 420];
exp_runlengths = [0.970, 1.310, 2.420, 1.660, 1.960, 2.86];
exp_err_runlengths = [0.180, 0.320, 0.350, 0.940, 0.310, 0.72];
exp_lifetimes = [1.8, 2.1, 7.1, 5.2, 8.3, 17.9];
exp_err_lifetimes = [0.6, 0.7, 1.7, 5.9, 2.6, 3.9];
exp_velocities = [600, 710, 360, 310, 310, 180];
exp_err_velocities = [75, 110, 50, 79, 40, 40];
%} 

color = [0, 0.4470, 0.7410; 0.9290, 0.6940, 0.1250; 0.2660, 0.7740, 0.3880];
exp_color = [0.85 0 0];

fig1 = figure();
set(fig1, 'Position', [50, 50, 1440, 2*285])
    
subplot(1, 4, 1)
hold all
exp_data = errorbar(exp_concs, exp_runlengths, exp_err_runlengths, '^', ...
    'MarkerSize', 12, 'LineWidth', 2, 'MarkerEdgeColor', exp_color);
exp_data.MarkerFaceColor = exp_data.MarkerEdgeColor;
exp_data.Color = exp_data.MarkerEdgeColor;
for i_run = 1 : n_runs
    sim_data = errorbar(concentrations, runlengths(i_run, :), err_runlengths(i_run, :), "o", ... 
         'MarkerSize', 7,'LineWidth', 2, 'MarkerEdgeColor', color(i_run, :));
    sim_data.MarkerFaceColor = sim_data.MarkerEdgeColor;
    sim_data.Color = sim_data.MarkerEdgeColor;
end
ylabel('Run length (microns)', 'FontSize', 14);
set(gca, 'FontSize', 14);
xlim([0 440]);
%ylim([0 4]);

subplot(1, 4, 2)
exp_data = errorbar(exp_concs, exp_lifetimes, exp_err_lifetimes, '^', ... 
    'MarkerSize', 12, 'LineWidth', 2, 'MarkerEdgeColor', exp_color);
exp_data.MarkerFaceColor = exp_data.MarkerEdgeColor;
exp_data.Color = exp_data.MarkerEdgeColor;
hold on
for i_run = 1 : n_runs
    sim_data = errorbar(concentrations, lifetimes(i_run, :), err_lifetimes(i_run, :), "o", ...
        'MarkerSize', 7, 'LineWidth', 2, 'MarkerEdgeColor', color(i_run, :));
    sim_data.MarkerFaceColor = sim_data.MarkerEdgeColor;
    sim_data.Color = sim_data.MarkerEdgeColor;
end
xlabel('KIF4A concentration (pM)', 'FontSize', 14);
ylabel('Life time (seconds)', 'FontSize', 14);
set(gca, 'FontSize', 14);
xlim([0 440]);
%ylim([-2 24]);

subplot(1, 4, 3)
exp_data = errorbar(exp_concs, exp_velocities, exp_err_velocities, '^', ...
    'MarkerSize', 12, 'LineWidth', 2, 'MarkerEdgeColor', exp_color);
exp_data.MarkerFaceColor = exp_data.MarkerEdgeColor;
exp_data.Color = exp_data.MarkerEdgeColor;
hold on
for i_run = 1 : n_runs
    sim_data = errorbar(concentrations, velocities(i_run, :), err_velocities(i_run, :), "o", ... 
       'MarkerSize', 7, 'LineWidth', 2, 'MarkerEdgeColor', color(i_run, :));
   sim_data.MarkerFaceColor = sim_data.MarkerEdgeColor;
   sim_data.Color = sim_data.MarkerEdgeColor;
end
ylabel('Velocity (nm/s)', 'FontSize', 14);
set(gca, 'FontSize', 14);
xlim([0 440]);
%ylim([0 900]);

subplot(1, 4, 4);
hold on
exp_data = errorbar([0, 0], [0, 0], '^', 'LineWidth', 2, 'MarkerSize', 12, ...
    'MarkerEdgeColor', exp_color);
exp_data.MarkerFaceColor = exp_data.MarkerEdgeColor;
exp_data.Color = exp_data.MarkerEdgeColor;
for i_run = 1 : n_runs
    sim_data = errorbar([0, 0], [0, 0], "o", 'LineWidth', 2, 'MarkerSize', 7, ...
        'MarkerEdgeColor', color(i_run, :));
    sim_data.MarkerFaceColor = sim_data.MarkerEdgeColor;
    sim_data.Color = sim_data.MarkerEdgeColor;
end
xlim([1 2]);
ylim([1 2]);
axis off
if length(baseNames) <= length(labelNames)
    leg = legend(["Experiment", labelNames], 'location', 'best', 'FontSize', 12);
else
    leg = legend(["Experiment", baseNames], 'location', 'northwest', 'FontSize', 12);
end
set(leg,'units','normalized');
set(leg,'position',[0.8,0.805,0.1,0.1])

