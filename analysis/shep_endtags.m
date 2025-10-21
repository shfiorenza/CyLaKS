%{
clear variables;
base_names = ["shep_0.1nM_1nM", "shep_0.1nM_10nM", "shep_0.1nM_100nM", ... 
    "shep_1nM_1000nM","shep_1nM_10nM", "shep_1nM_50nM", "shep_1nM_100nM", ...
    "shep_1nM_250nM", "shep_1nM_500nM", "shep_1nM_1000nM"];
%base_names = ["shep_1nM_100nM"];
folder = "out_final";

mt_lengths = [500, 750, 1000, 1250]; % in n_sites
seeds = [0, 1, 2, 3, 4, 5];

dir = sprintf("../%s", folder);

site_size = 0.0082; % in um
n_runs = length(base_names);
n_mts = length(mt_lengths);
n_seeds = length(seeds);

endtag_lengths = zeros(n_runs, n_mts, n_seeds);
avg_endtag_length = zeros(n_runs, n_mts);
err_endtag_length = zeros(n_runs, n_mts);
for i_run = 1 : n_runs
    for i_mt = 1 : n_mts
        for i_seed = 1 : n_seeds
            sim_name = sprintf("%s/%s_8_%i_0.6kT_3x_5x_%i", dir, base_names(i_run), mt_lengths(i_mt), seeds(i_seed))
            endtag_lengths(i_run, i_mt, i_seed) = get_endtag_length_xlink(sim_name);
        end
    end
    for i_mt = 1 : n_mts
        avg_endtag_length(i_run, i_mt) = mean(endtag_lengths(i_run, i_mt, :));
        var = 0.0;
        for i_seed = 1 : n_seeds
            var = var + double(avg_endtag_length(i_run, i_mt) - endtag_lengths(i_run, i_mt, i_seed))^2 / (n_seeds - 1);
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

fig1 = figure();
set(fig1, 'Position', [50, 50, 960, 540])
hold all;
% Plot sim data
%{
for i_run = 1 : n_runs
    sim_data = errorbar(mt_lengths * site_size, avg_endtag_length(i_run, :, :), ...
        err_endtag_length(i_run, :), 'o','MarkerSize', 12, 'LineWidth', 2);%, ...
        %'MarkerEdgeColor', color(i_run, :));
   %sim_data.MarkerFaceColor = sim_data.MarkerEdgeColor;
   %sim_data.Color = sim_data.MarkerFaceColor;
end
%}
markers = ['^', '^', '^', 'o', 'o', 'o', 'o', 'o', 'o','diamond'];
colors = [255 160 160; 247 78 214; 159 4 77; 87 87 249; 0 0 192; 65 240 174; 0 128 0; 192 192 0; 255 128 0; 249 64 64]/255;
for i_seed = 1 : n_seeds
    for i_run = 1 : n_runs
        data = plot(mt_lengths * site_size, endtag_lengths(i_run, :, i_seed), markers(i_run),'MarkerSize', 8, 'LineWidth', 2, 'MarkerEdgeColor', colors(i_run, :));
        data.MarkerFaceColor = data.MarkerEdgeColor;
        %data.Color = data.MarkerEdgeColor;
    end
end
%}

%{
for i_mt = 1 : n_mts
    sim_data = plot(mt_lengths * site_size, avg_endtag_length(:, i_mt), ...
        'o','MarkerSize', 14, 'MarkerEdgeColor', color(i_mt, :));
    sim_data.MarkerFaceColor = sim_data.MarkerEdgeColor;
    sim_data.Color = sim_data.MarkerFaceColor;
end
%}

%xlabel("Microtubule length (\mum)", 'FontSize', 18);
xlabel("Range of potential (\mum)");
ylabel("Endtag length (\mum)");
set(gca, 'FontSize', 22);
%ylim([0 1.9]); % 2]); %-0.25 11]);
%yticks([0 0.5 1 1.5 2]); %5 10]);
xlim([0 12]);
%ylim([0 6]);
%xticks([0.1 1 10]);
%xticklabels([0.1 1 10]); 
%set(gca, 'XScale', 'log');


% Varied conc stylistic stuff

xlabel("Microtubule length (\mum)", 'FontSize', 18);
ylabel("End-zone length (\mum)", 'Fontsize', 18);
legendLabel = ["0.1, 1", "0.1, 10", "0.1, 100", ...
               "1, 1", "1, 10", "1, 50", "1, 100", "1, 250", "1, 500", "1, 1000"];
leg = legend(legendLabel,'location', 'northeastoutside', 'FontSize', 22);
legend('boxoff');
title(leg, 'PRC1, K401 (nM)')
ylim([0 6]);
xlim([2 11]);
set(gca,'box','off')
set(gca, 'FontSize', 24);
set(gca,'TickDir','out');
set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);