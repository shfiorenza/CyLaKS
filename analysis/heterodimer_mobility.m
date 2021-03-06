%{
clear variables;
baseNames = ["hetero_tubulin"];
folder = "run_tubulin_mobility";

dir = sprintf("/home/shane/projects/CyLaKS/%s", folder);

color_exp = [0, 0.4470, 0.7410]; %[1 0 0]; %[0.55 0.55 0.55];

p_hetero = [0.0, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0];
seeds = [0, 1, 2, 3, 4, 5];

n_fractions = length(p_hetero);

runlengths = zeros(n_fractions, 1);
err_runlengths = zeros(n_fractions, 1);
lifetimes = zeros(n_fractions, 1);
err_lifetimes = zeros(n_fractions, 1);
velocities = zeros(n_fractions, 1);
err_velocities = zeros(n_fractions, 1);

for i_frac = 1 : n_fractions
    simName = sprintf("%s/%s_%g", dir, baseNames(1), p_hetero(i_frac));
    if length(seeds) > 1
        mot_stats = get_motor_stats(simName, seeds);
    else
        mot_stats = get_motor_stats(simName);
    end
    runlengths(i_frac) = mot_stats(1);
    err_runlengths(i_frac) = mot_stats(2);
    lifetimes(i_frac) = mot_stats(3);
    err_lifetimes(i_frac) = mot_stats(4);
    velocities(i_frac) = mot_stats(5);
    err_velocities(i_frac) = mot_stats(6);
end
%}
fig1 = figure();
set(fig1, 'Position', [50, 50, 1400, 350])  % For 3 subplots

subplot(1, 3, 1) 
hold on
set(gca, 'FontSize', 20);
exp_data = errorbar(p_hetero, runlengths/1000, err_runlengths/1000, 'o', ...
    'MarkerSize', 8, 'LineWidth', 2, 'MarkerEdgeColor', color_exp);
exp_data.MarkerFaceColor = exp_data.MarkerEdgeColor;
exp_data.Color = exp_data.MarkerEdgeColor;
ylabel('Run Length (\mum)', 'FontSize', 20);
xlim([-0.05 1.05]);
ylim([0 13]);
yticks([0 4 8 12]);

subplot(1, 3, 2 ) %4, 2)
hold on
set(gca, 'FontSize', 20);
exp_data = errorbar(p_hetero, lifetimes, err_lifetimes, 'o', ... 
    'MarkerSize', 8, 'LineWidth', 2, 'MarkerEdgeColor', color_exp);
exp_data.MarkerFaceColor = exp_data.MarkerEdgeColor;
exp_data.Color = exp_data.MarkerEdgeColor;
xlabel('Fraction of sites that are mutant', 'FontSize', 20);
ylabel('Lifetime (s)', 'FontSize', 20);
xlim([-0.05 1.05]);
ylim([0 23]);
yticks([0 7 14 21]);

subplot(1, 3, 3) %4, 3)
hold on
set(gca, 'FontSize', 20);
exp_data = errorbar(p_hetero, velocities/1000, err_velocities/1000, 'o', ...
    'MarkerSize', 8, 'LineWidth', 2, 'MarkerEdgeColor', color_exp);
exp_data.MarkerFaceColor = exp_data.MarkerEdgeColor;
exp_data.Color = exp_data.MarkerEdgeColor;
ylabel('Velocity (\mum/s)', 'FontSize', 20);
xlim([-0.05 1.05]);
ylim([0 1]);
yticks([0 0.3 0.6 0.9]);

%suptitle('Kinesin mobility with mutant sites (k_D = 4k^0_D) randomly inserted')