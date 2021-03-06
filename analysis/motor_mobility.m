%
clear variables;
baseNames = ["test"];
folder = "";
% Data for baseline mobility w/o any coop -- Kif4A
%{
baseNames = ["mobility_baseline_Kif4A"];
folder = "run_mobility_baseline_Kif4A";
%}
% Data for baseline mobility w/o any coop -- k401 (use B set of exp data)
%{
baseNames = ["mobility_baseline_K401"];
folder = "run_mobility_baseline_K401";
%}
% Data for short-range coop
%{
baseNames = ["mobility_short_8"];
folder = "run_mobility_short";
%}
% Data for short- & long-range coop that doesn't affect stepping
%{
baseNames = ["mobility_both_nostep"];
folder = "run_mobility_both_nostep";
%}
% Data for short- & long-range coop that also affects stepping
%{
baseNames = ["mobility_both"];
folder = "run_mobility_both";
%}
% Data for all short-range coop runs
%{
baseNames = [];
for i_energy = 0 : 2 : 10
   baseNames = [baseNames sprintf("mobility_short_%i", i_energy)]; 
end
folder = "run_mobility_short";
%}
% Data for long-range coop ONLY 
%{
baseNames = ["mobility_long", "../run_mobility_both/mobility_both"];
folder = "run_mobility_long";
%}

concentrations = [20, 50, 80, 120, 220, 420];
concs_index = [1, 2, 3, 4, 5, 6];
seeds = [0]; % , 1, 2, 3];

dir = sprintf("/home/shane/projects/LAT-CAT/%s", folder);

% Kif4A data
exp_runlengths = [970, 1310, 2420, 1660, 1960, 2860];
exp_err_runlengths = [180, 320, 350, 940, 310, 720];
exp_lifetimes = [1.8, 2.1, 7.1, 5.2, 8.3, 17.9];
exp_err_lifetimes = [0.6, 0.7, 1.7, 5.9, 2.6, 3.9];
exp_velocities = [600, 710, 360, 310, 310, 180];
exp_err_velocities = [75, 110, 50, 79, 40, 40];
% K401 data
exp_runlengthsB = [1058, 1255, 1612, 2720, 1945, 3471];
exp_err_runlengthsB = [124, 99, 160, 652, 283, 480];
exp_lifetimesB = [1.944, 2.684, 5.724, 7.864, 5.303, 15.69];
exp_err_lifetimesB = [0.35, 0.4, 0.95, 0.55, 0.65, 2.2];
exp_velocitiesB = [659.3, 590.9, 270.8, 376.5, 345.9, 210.3];
exp_err_velocitiesB = [25.2, 18.9, 23, 27.3, 21, 31.5];

n_runs = length(baseNames);
n_concs = length(concentrations);
runlengths = zeros(n_runs, n_concs);
err_runlengths = zeros(n_runs, n_concs);
lifetimes = zeros(n_runs, n_concs);
err_lifetimes = zeros(n_runs, n_concs);
velocities = zeros(n_runs, n_concs);
err_velocities = zeros(n_runs, n_concs);
for i_run = 1 : n_runs
    for i_concs = 1 : n_concs
        conc = int32(concentrations(i_concs));
        simName = sprintf("%s/%s_%i", dir, baseNames(i_run), conc);
        if length(seeds) > 1
            mot_stats = get_motor_stats(simName, seeds);
        else
            mot_stats = get_motor_stats(simName);
        end
        runlengths(i_run, i_concs) = mot_stats(1);
        err_runlengths(i_run, i_concs) = mot_stats(2);
        lifetimes(i_run, i_concs) = mot_stats(3);
        err_lifetimes(i_run, i_concs) = mot_stats(4);
        velocities(i_run, i_concs) = mot_stats(5);
        err_velocities(i_run, i_concs) = mot_stats(6);
    end
end
%}

%color = [0 0.447 0.741; 0.85, 0.325, 0.098; 0.929, 0.694, 0.125; ...
%    0.494, 0.184, 0.556; 0.466, 0.674, 0.188; 0.301, 0.745, 0.933];
color = [0 0.447 0.741; 0.266, 0.674, 0.388];
color_base = [0, 0.4470, 0.7410];
color_exp = [1 0 0]; %[0.55 0.55 0.55];
color_run = [51 0 255] / 255;
color_life = [255 153 51] / 255;
color_vel = [102 0 104] / 255;

i_start = 1;
%n_runs = 3;

fig1 = figure();
set(fig1, 'Position', [50, 50, 1.25*1040, 1.25*285])  % For 3 subplots
%set(fig1, 'Position', [50, 50, 1440, 285])  % For 4 subplots

subplot(1, 3, 1) %4, 1)
hold all;
exp_data = errorbar(concs_index, exp_runlengths, exp_err_runlengths, '^', ...
    'MarkerSize', 16, 'LineWidth', 2, 'MarkerEdgeColor', color_exp);
exp_data.MarkerFaceColor = exp_data.MarkerEdgeColor;
exp_data.Color = exp_data.MarkerEdgeColor;
for i_run = i_start : n_runs
    sim_data = errorbar(concs_index, runlengths(i_run, :), err_runlengths(i_run, :), "o", ... 
         'MarkerSize', 12,'LineWidth', 2, 'MarkerEdgeColor', color(i_run, :)); %color_run);
    sim_data.MarkerFaceColor = sim_data.MarkerEdgeColor;
    sim_data.Color = sim_data.MarkerEdgeColor;
end
ylabel('Run Length (nm)', 'FontSize', 14);
set(gca, 'FontSize', 14);
xlim([0 7]);
xticks(1 : 6);
xticklabels({'20', '50', '80', '120', '220', '420'});
%xticklabels({''});
xtickangle(45);
ylim([0 4000]); % 5000]); % for B-set of exp data
yticks([0 1500 3000]); % 2000 4000]); % for B set of exp data

legendLabel = ["Experiment", "Simulation"];
leg = legend(legendLabel,'location', 'northwest', 'FontSize', 12);
legend('boxoff');
set(leg,'units','normalized');
set(leg,'position',[0.14,0.8,0.1,0.1])
%}

subplot(1, 3, 2 ) %4, 2)
hold all;
exp_data = errorbar(concs_index, exp_lifetimes, exp_err_lifetimes, '^', ... 
    'MarkerSize', 16, 'LineWidth', 2, 'MarkerEdgeColor', color_exp);
exp_data.MarkerFaceColor = exp_data.MarkerEdgeColor;
exp_data.Color = exp_data.MarkerEdgeColor;
for i_run = i_start : n_runs
    sim_data = errorbar(concs_index, lifetimes(i_run, :), err_lifetimes(i_run, :), "o", ...
        'MarkerSize', 12, 'LineWidth', 2, 'MarkerEdgeColor', color(i_run, :)); %color_life);
    sim_data.MarkerFaceColor = sim_data.MarkerEdgeColor;
    sim_data.Color = sim_data.MarkerEdgeColor;
end
xlabel('Kif4A Concentration (pM)', 'FontSize', 14);
ylabel('Lifetime (s)', 'FontSize', 14);
set(gca, 'FontSize', 14);
xlim([0 7]);
xticks(1 : 6);
xticklabels({'20', '50', '80', '120', '220', '420'});
%xticklabels({''});
xtickangle(45);
ylim([-2 24]);
yticks([0 10 20]);

subplot(1, 3, 3) %4, 3)
hold all;
exp_data = errorbar(concs_index, exp_velocities, exp_err_velocities, '^', ...
    'MarkerSize', 16, 'LineWidth', 2, 'MarkerEdgeColor', color_exp);
exp_data.MarkerFaceColor = exp_data.MarkerEdgeColor;
exp_data.Color = exp_data.MarkerEdgeColor;
for i_run = i_start : n_runs
    sim_data = errorbar(concs_index, velocities(i_run, :), err_velocities(i_run, :), "o", ... 
       'MarkerSize', 12, 'LineWidth', 2, 'MarkerEdgeColor', color(i_run, :)); %color_vel);
   sim_data.MarkerFaceColor = sim_data.MarkerEdgeColor;
   sim_data.Color = sim_data.MarkerEdgeColor;
end
ylabel('Velocity (nm/s)', 'FontSize', 14);
set(gca, 'FontSize', 14);
xlim([0 7]);
xticks(1 : 6);
xticklabels({'20', '50', '80', '120', '220', '420'});
%xticklabels({''});
xtickangle(45);
ylim([0 900]);
yticks([0 400 800]);

% Add a fourth 'ghost plot' for legends when needed 
%{
subplot(1, 4, 4);
hold all;
exp_data = errorbar([0, 0], [0, 0], '^', 'LineWidth', 2, 'MarkerSize', 16, ...
    'MarkerEdgeColor', color_exp);
exp_data.MarkerFaceColor = exp_data.MarkerEdgeColor;
exp_data.Color = exp_data.MarkerEdgeColor;
for i_run = i_start : n_runs
    sim_data = errorbar([0, 0], [0, 0], "o", 'LineWidth', 2, 'MarkerSize', 12, ...
        'MarkerEdgeColor', color(i_run, :));
    sim_data.MarkerFaceColor = sim_data.MarkerEdgeColor;
    sim_data.Color = sim_data.MarkerEdgeColor;
end
xlim([1 2]);
ylim([1 2]);
axis off
% Style stuff for shortScan data -- ensure subplot() calls above change
%{
legendLabel = ["Experiment", "E = 0 KT"];
for i_energy = 2 : 2 : 10
    legendLabel = [legendLabel sprintf("E = -%i kT", i_energy)];
end
%}
% Style stuff for long_only vs both stuff
legendLabel = ["Experiment", "Simulation"];
leg = legend(legendLabel,  'location', 'northwest', 'FontSize', 12);
set(leg, 'units', 'normalized');
set(leg, 'position',[0.735,0.775,0.1,0.1]);
legend('boxoff');
%}
