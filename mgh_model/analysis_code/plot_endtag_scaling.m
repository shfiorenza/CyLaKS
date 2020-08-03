%{
clear variables;

base_names = ["data_endtag/endtag", "data_endtag_short/endtag_shortCoopB"]; %, "data_endtag_fullCoop/endtag"];
mt_lengths = [250, 500, 750, 1000, 1250, 1750]; % in n_sites
seeds = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];

dir = "/home/shane/Projects/overlap_analysis/mgh_model";

exp_mt_lengths = [2.39, 3.99, 5.57, 7.15, 8.75, 10.33, 13.51]; % in um
exp_err_mt_lengths =  [0.79, 0.79, 0.79, 0.79, 0.79, 0.79, 0.79]; % in um;
exp_endtag_lengths = [1.19, 1.27, 1.43, 1.61, 1.63, 1.74, 2.21]; % in um
exp_err_endtag_lengths = [0.11, 0.10, 0.12, 0.12, 0.18, 0.16, 0.18, ]; % in um

site_size = 0.008; % in um
n_runs = length(base_names);
n_mts = length(mt_lengths);
n_seeds = length(seeds);

avg_endtag_length = zeros(n_runs, n_mts);
err_endtag_length = zeros(n_runs, n_mts);
for i_run = 1 : n_runs
    endtag_lengths = zeros(n_mts, n_seeds);
    for i_mt = 1 : n_mts
        for i_seed = 1 : n_seeds
            sim_name = sprintf("%s/%s_%i_%i", dir, base_names(i_run), mt_lengths(i_mt), seeds(i_seed));
            endtag_lengths(i_mt, i_seed) = get_endtag_length(sim_name);
        end
    end
    for i_mt = 1 : n_mts
        avg_endtag_length(i_run, i_mt) = mean(endtag_lengths(i_mt, :));
        var = 0.0;
        for i_seed = 1 : n_seeds
            var = var + double(avg_endtag_length(i_mt) - endtag_lengths(i_mt, i_seed))^2 / (n_seeds - 1);
        end
        err_endtag_length(i_run, i_mt) = sqrt(var / n_seeds);
    end
end

%}

%colors = ["b", "g"];
color = [0.929, 0.694, 0.125; 0.494, 0.184, 0.556; 0.0, 0.4470, 0.7410]; 

fig1 = figure();
set(fig1, 'Position', [50, 50, 720, 720])
hold all;

% Plot exp data with just vertical error bars
exp_data = errorbar(exp_mt_lengths, exp_endtag_lengths, ...
    exp_err_endtag_lengths, 'r^', 'MarkerSize', 12);
exp_data.MarkerFaceColor = exp_data.MarkerEdgeColor;
%exp_data.Color = exp_data.MarkerEdgeColor;
% Plot sim data
for i_run = 1 : n_runs
    sim_data = errorbar(mt_lengths * site_size, avg_endtag_length(i_run, :), ...
        err_endtag_length(i_run, :),'o','MarkerSize', 12, 'LineWidth', 2, ...
        'MarkerEdgeColor', color(i_run, :));
    sim_data.MarkerFaceColor = sim_data.MarkerEdgeColor;
    sim_data.Color = sim_data.MarkerEdgeColor;
end
% Use 3rd party errorbarxy.m to plot horizontal error bars for exp data 
errorbarxy(exp_mt_lengths, exp_endtag_lengths, exp_err_mt_lengths, ...
    exp_err_endtag_lengths, {'r^', 'r', 'r'});

xlabel("Microtubule length (microns)", 'FontSize', 18);
ylabel("Endtag length (microns)", 'Fontsize', 18);
set(gca, 'FontSize', 18);
legend(["Experiment", "Sim", "Sim - short-range coop"], ... , "Sim - long-range coop"], ...
    'location', 'northwest', 'FontSize', 18);
%ylim([0 3]);
xlim([0 16]);
