%{
clear variables;
base_names = ["endtag"];
folder = "run_endtag_vs_coop";
% Data for initial plot that compares to previous paper
%{
base_names = ["endtag_1.5nM"];
folder = "run_endtag_baseline";
%}
% Data for varied conc runs 
%{
base_names = ["endtag_0.02nM", "endtag_1nM", "endtag_2nM", "endtag_4nM", "endtag_6nM"];
folder = "run_endtag_baseline";
%}
% Data for 'final' results using both short-range & long-range that affects stepping
%{
base_names = ["endtag_both"];
folder = "run_endtag_both";
%}
% Data for 'final' varied conc results
%{
base_names = ["endtag_final_0.02nM", "endtag_final_1nM", "endtag_final_2nM", ...
    "endtag_final_4nM", "endtag_final_6nM"];
folder = "run_endtag_final";
%}
% Data for all short-range coop runs
%{
base_names = [];
for i_energy = 0 : 2 : 10
   base_names = [base_names sprintf("endtag_short_%i", i_energy)]; 
end
folder = "run_endtag_short";
%}
% Data for 10x processivity w/ no coop runs
%{
base_names = ["endtag_10xProc"];
folder = "run_endtag_10xProc"; 
%}
% Data for end-tags with short- & long-range coop but no stepping FX
%{
base_names = ["endtag_both_nostep", "../run_endtag_both/endtag_both"];
folder = "run_endtag_both_nostep";
%}
% Data for full model, sans stepping FX, sans long range, sans short range
%{
base_names = ["endtag_both", "../run_endtag_both_nostep/endtag_both_nostep", ...
    "../run_endtag_both_nobind/endtag_both_nobind", "../run_endtag_long/endtag_long"]; 
folder = "run_endtag_both";
%}

mt_lengths = [250, 500, 750, 1000, 1250, 1750]; % in n_sites
ranges = [10, 50, 100, 1000];
seeds = [0, 1, 2, 3]; %, 4, 5, 6, 7, 8, 9];

dir = sprintf("/home/shane/projects/CyLaKS/%s", folder);

exp_mt_lengths = [2.39, 3.99, 5.57, 7.15, 8.75, 10.33, 13.51]; % in um
exp_err_mt_lengths =  [0.79, 0.79, 0.79, 0.79, 0.79, 0.79, 0.79]; % in um;
exp_endtag_lengths = [1.19, 1.27, 1.43, 1.61, 1.63, 1.74, 2.21]; % in um
exp_err_endtag_lengths = [0.11, 0.10, 0.12, 0.12, 0.18, 0.16, 0.18, ]; % in um

site_size = 0.008; % in um
%n_runs = length(base_names);
n_runs = length(ranges);
n_mts = length(mt_lengths);
n_seeds = length(seeds);

avg_endtag_length = zeros(n_runs, n_mts);
err_endtag_length = zeros(n_runs, n_mts);
for i_run = 1 : n_runs
    endtag_lengths = zeros(n_mts, n_seeds);
    for i_mt = 1 : n_mts
        for i_seed = 1 : n_seeds
            %sim_name = sprintf("%s/%s_%i_%i", dir, base_names(i_run),conc, seeds(i_seed));
            %sim_name = sprintf("%s/%s_%i_%i", dir, base_names(i_run), mt_lengths(i_mt), seeds(i_seed));
            %for range = i_range : length(ranges)
            sim_name = sprintf("%s/%s_%i_%i_%i", dir, base_names(1), mt_lengths(i_mt),  ...
                ranges(i_run), seeds(i_seed))
           % end
            endtag_lengths(i_mt, i_seed) = get_endtag_length(sim_name);
        end
    end
    for i_mt = 1 : n_mts
        avg_endtag_length(i_run, i_mt) = mean(endtag_lengths(i_mt, :));
        var = 0.0;
        for i_seed = 1 : n_seeds
            var = var + double(avg_endtag_length(i_run, i_mt) - endtag_lengths(i_mt, i_seed))^2 / (n_seeds - 1);
        end
        err_endtag_length(i_run, i_mt) = sqrt(var / n_seeds);
    end
end
%}

color = [0 0.447 0.741; 0.85, 0.325, 0.098; 0.929, 0.694, 0.125; ...
    0.494, 0.184, 0.556; 0.466, 0.674, 0.188; 0.301, 0.745, 0.933];
marker = {'o', 'o', 'o', 'o', 'o'};
% Color & markers for baseline varied conc
%{
color = [160 160 160; 128 128 128; 96 96 80; 64 64 64; 0 0 0 ] / 255;
marker = ['o', 's', '^', 'v', 'd'];
%}

%fig1 = figure();
%set(fig1, 'Position', [50, 50, 1080, 720])
hold all;
%{
% Plot exp data with just vertical error bars
exp_data = errorbar(exp_mt_lengths, exp_endtag_lengths, ...
    exp_err_endtag_lengths, 'r^', 'MarkerSize', 16);
exp_data.MarkerFaceColor = exp_data.MarkerEdgeColor;
% Plot sim data
%}

for i_run = 1 : n_runs
    sim_data = errorbar(mt_lengths * site_size, avg_endtag_length(i_run, :), ...
        err_endtag_length(i_run, :), 'o','MarkerSize', 12, 'LineWidth', 2, ...
        'MarkerEdgeColor', color(i_run, :));
   sim_data.MarkerFaceColor = sim_data.MarkerEdgeColor;
   sim_data.Color = sim_data.MarkerFaceColor;
end
%}
%{
for i_mt = 1 : n_mts
    sim_data = errorbar(ranges * site_size, avg_endtag_length(:, i_mt), ...
        err_endtag_length(:, i_mt), 'o','MarkerSize', 12, 'LineWidth', 2, ...
        'MarkerEdgeColor', color(i_mt, :));
   sim_data.MarkerFaceColor = sim_data.MarkerEdgeColor;
   sim_data.Color = sim_data.MarkerFaceColor;
end
%}
%{
% Use 3rd party errorbarxy.m to plot horizontal error bars for exp data 
errorbarxy(exp_mt_lengths, exp_endtag_lengths, exp_err_mt_lengths, ...
    exp_err_endtag_lengths, {'r^', 'r', 'r'});
%}

% Plot for varied conc ratio
%{
kif4a_concs = [0.002, 1, 2, 4, 6];
n_concs = length(kif4a_concs);
% Calculate average endtag-to-microtubule length ratio for each conc
avg_length_ratio = zeros(n_concs, 1);
err_length_ratio = zeros(n_concs, 1);
for i_conc = 1 : n_concs
    for i_mt = 1 : n_mts
       ratio = avg_endtag_length(i_conc, i_mt) / (mt_lengths(i_mt) * site_size);
       avg_length_ratio(i_conc, 1) = avg_length_ratio(i_conc, 1) + ratio / n_mts; 
    end
    % calculate error of average
    var = 0.0;
    for i_mt = 1 : n_mts
        ratio = avg_endtag_length(i_conc, i_mt) / (mt_lengths(i_mt) * site_size);
        delta = avg_length_ratio(i_conc, 1) - ratio;
        var = var + (delta * delta / (n_mts - 1));
    end
    err_length_ratio(i_conc, 1) = sqrt(var / n_mts);
end
% Plot that jawn
for i_conc = 1 : n_concs
    sim_data = errorbar(kif4a_concs(i_conc), avg_length_ratio(i_conc), err_length_ratio(i_conc), ...
        marker(i_conc),'MarkerSize', 12, 'LineWidth', 2, 'MarkerEdgeColor', color(i_conc, :));
    sim_data.MarkerFaceColor = sim_data.MarkerEdgeColor;
    sim_data.Color = sim_data.MarkerFaceColor;
end
%}


xlabel("Microtubule length (\mum)", 'FontSize', 18);
ylabel("Endtag length (\mum)", 'Fontsize', 18);
set(gca, 'FontSize', 18);
%legendLabel = ["Experiment", "Simulation"];
%legend(legendLabel,'location', 'northwest', 'FontSize', 18);
%legend('boxoff');
ylim([0 3.5]); % 2]); %-0.25 11]);
yticks([0 1 2 3]); %5 10]);
xlim([0 15]); % 12]);
xticks([0 5 10 15]);

% Varied conc stylistic stuff
%{
xlabel("Microtubule length, L_{M} (\mum)", 'FontSize', 18);
ylabel("Endtag length, L_{ET} (\mum)", 'Fontsize', 18);
legendLabel = ["0.02 nM", "1", "2", "4", "6"];
legend(legendLabel,'location', 'northwest', 'FontSize', 16);
legend('boxoff');
%}

% Varied conc ratio stylistic stuff
%{
xlabel("Kif4A concentration (nM)", 'FontSize', 18);
ylabel("L_{ET} / L_{M}", 'Fontsize', 18);
set(gca, 'FontSize', 18);
ylim([0.0 0.8]);
yticks([0.0 0.4 0.8]);
xlim([0 6.5]);
xticks([0 2 4 6]);
%}

% Short-range coop stylistic stuff
%{
legendLabel = ["Experiment", "E = 0 KT"];
for i_energy = 2 : 2 : 10
    legendLabel = [legendLabel sprintf("E = -%i kT", i_energy)];
end
legend(legendLabel,'location', 'northwest', 'FontSize', 16);
legend('boxoff');
%}

% 10x processivity w/o cooperativity stylistic stuff
%{
legendLabel = ["Experiment", "Simulation w/ 10x processivity"];
legend(legendLabel,'location', 'northwest', 'FontSize', 18);
legend('boxoff');
%}

% Short- & long-range coop w/o stepping FX stylistic stuff
legendLabel = ["Range = 0.08 \mum", "Range = 0.4 \mum", "Range = 0.8 \mum", "Range = 8 \mum"];
legend(legendLabel,'location', 'northwest', 'FontSize', 18);
legend('boxoff');
%{
legendLabel = ["Experiment", "Simulation", ...
    "Simulation - no long-range stepping cooperativity", ...
    "Simulation - no long-range binding cooperativity", ...
    "Simulation - no short-range binding cooperativity"];
legend(legendLabel,'location', 'northwest', 'FontSize', 18);
legend('boxoff');
%}
