concs_index = [1, 2, 3, 4, 5, 6, 7, 8];
color_exp = [0, 0.4470, 0.7410]; %[1 0 0]; %[0.55 0.55 0.55];
exp_runlengths = [1.3, 1.0 ,0.9, 0.6, 0.7, 0.4, 0.2, 0.1];
exp_err_runlengths = [0.5, 0.4, 0.4, 0.2, 0.3, 0.2, 0.09, 0.1];
exp_lifetimes = [2.4, 1.6, 2.0, 1.6, 1.2, 0.7, 0.3, 0.4];
exp_err_lifetimes = [0.9, 0.6, 0.8, 0.3, 0.5, 0.3, 0.1, 0.2];
exp_velocities = [570, 600, 560, 600, 580, 600, 640, 570];
exp_err_velocities = [60, 20,60, 40, 40, 50, 80, 40];

fig1 = figure();
set(fig1, 'Position', [50, 50, 1000, 700])  % For 3 subplots

subplot(1, 3, 1) 
exp_data = errorbar(concs_index, exp_runlengths, exp_err_runlengths, '^', ...
    'MarkerSize', 16, 'LineWidth', 2, 'MarkerEdgeColor', color_exp);
exp_data.MarkerFaceColor = exp_data.MarkerEdgeColor;
exp_data.Color = exp_data.MarkerEdgeColor;
ylabel('Run Length (nm)', 'FontSize', 14);
set(gca, 'FontSize', 14);
xlim([0 length(concs_index)+1]);
xticks(1 : length(concs_index));
xticklabels({'0', '5', '10', '20', '40', '60', '80', '100'});
xtickangle(45);
ylim([0 2]); % 5000]); % for B-set of exp data
yticks([0 0.5 1 1.5 2]); % 2000 4000]); % for B set of exp data


subplot(1, 3, 2 ) %4, 2)
exp_data = errorbar(concs_index, exp_lifetimes, exp_err_lifetimes, '^', ... 
    'MarkerSize', 16, 'LineWidth', 2, 'MarkerEdgeColor', color_exp);
exp_data.MarkerFaceColor = exp_data.MarkerEdgeColor;
exp_data.Color = exp_data.MarkerEdgeColor;
xlabel('Mutant tubulin percentage', 'FontSize', 14);
ylabel('Lifetime (s)', 'FontSize', 14);
set(gca, 'FontSize', 14);
xlim([0 length(concs_index)+1]);
xticks(1 : length(concs_index));
xticklabels({'0', '5', '10', '20', '40', '60', '80', '100'});
xtickangle(45);
ylim([0 4]);
yticks([0 1 2 3 4]);

subplot(1, 3, 3) %4, 3)
exp_data = errorbar(concs_index, exp_velocities, exp_err_velocities, '^', ...
    'MarkerSize', 16, 'LineWidth', 2, 'MarkerEdgeColor', color_exp);
exp_data.MarkerFaceColor = exp_data.MarkerEdgeColor;
exp_data.Color = exp_data.MarkerEdgeColor;
ylabel('Velocity (nm/s)', 'FontSize', 14);
set(gca, 'FontSize', 14);
xlim([0 length(concs_index)+1]);
xticks(1 : length(concs_index));
xticklabels({'0', '5', '10', '20', '40', '60', '80', '100'});
xtickangle(45);
ylim([0 800]);
yticks([0 200 400 600 800]);

suptitle('Kinesin mobility with mutant sites (k_D = 4k^0_D) randomly inserted')