%{
clear variables;
baseNames = ["hetero_tubulin"];
folder = "run_hetero_tubulin";

dir = sprintf("/home/shane/projects/CyLaKS/%s", folder);

color_exp = [0, 0.4470, 0.7410]; %[1 0 0]; %[0.55 0.55 0.55];

p_hetero = [0.0, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95, 1.0];
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
exp_data = errorbar(p_hetero, runlengths, err_runlengths, 'o', ...
    'MarkerSize', 8, 'LineWidth', 2, 'MarkerEdgeColor', color_exp);
exp_data.MarkerFaceColor = exp_data.MarkerEdgeColor;
exp_data.Color = exp_data.MarkerEdgeColor;
ylabel('Run Length (nm)', 'FontSize', 14);
set(gca, 'FontSize', 14);
xlim([-0.05 1.05]);
ylim([0 1800]);

subplot(1, 3, 2 ) %4, 2)
exp_data = errorbar(p_hetero, lifetimes, err_lifetimes, 'o', ... 
    'MarkerSize', 8, 'LineWidth', 2, 'MarkerEdgeColor', color_exp);
exp_data.MarkerFaceColor = exp_data.MarkerEdgeColor;
exp_data.Color = exp_data.MarkerEdgeColor;
xlabel('Fraction of sites that are mutant', 'FontSize', 14);
ylabel('Lifetime (s)', 'FontSize', 14);
set(gca, 'FontSize', 14);
xlim([-0.05 1.05]);
ylim([0 3]);

subplot(1, 3, 3) %4, 3)
exp_data = errorbar(p_hetero, velocities, err_velocities, 'o', ...
    'MarkerSize', 8, 'LineWidth', 2, 'MarkerEdgeColor', color_exp);
exp_data.MarkerFaceColor = exp_data.MarkerEdgeColor;
exp_data.Color = exp_data.MarkerEdgeColor;
ylabel('Velocity (nm/s)', 'FontSize', 14);
set(gca, 'FontSize', 14);
xlim([-0.05 1.05]);
ylim([0 700]);

%suptitle('Kinesin mobility with mutant sites (k_D = 4k^0_D) randomly inserted')