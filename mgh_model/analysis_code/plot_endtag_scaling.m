
clear variables;

base_name = "endtag";
mt_lengths = [250, 500, 750, 1000, 1250, 1750]; % in n_sites
seeds = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];

n_mts = length(mt_lengths);
n_seeds = length(seeds);
site_size = 0.008; % in um

exp_mt_lengths = [2.39, 3.99, 5.57, 7.15, 8.75, 10.33, 13.51]; % in um
exp_err_mt_lengths =  [0.79, 0.79, 0.79, 0.79, 0.79, 0.79, 0.79]; % in um;
exp_endtag_lengths = [1.19, 1.27, 1.43, 1.61, 1.63, 1.74, 2.21]; % in um
exp_err_endtag_lengths = [0.11, 0.10, 0.12, 0.12, 0.18, 0.16, 0.18, ]; % in um

dir = "/home/shane/Projects/overlap_analysis/mgh_model/data";

endtag_lengths = zeros(n_mts, n_seeds);
for i_mt = 1 : n_mts
   for i_seed = 1 : n_seeds
      sim_name = sprintf("%s/%s_%i_%i", dir, base_name, mt_lengths(i_mt), seeds(i_seed));
      endtag_lengths(i_mt, i_seed) = get_endtag_length(sim_name);
   end
end

avg_endtag_length = zeros(n_mts, 1);
err_endtag_length = zeros(n_mts, 1);
for i_mt = 1 : n_mts
   avg_endtag_length(i_mt) = mean(endtag_lengths(i_mt, :));
   var = 0.0;
   for i_seed = 1 : n_seeds
       var = var + double(avg_endtag_length(i_mt) - endtag_lengths(i_mt, i_seed))^2 / (n_seeds - 1);
   end
   err_endtag_length(i_mt) = sqrt(var / n_seeds);
end
%}

fig1 = figure();
set(fig1, 'Position', [50, 50, 480, 480])
hold all;
% Plot sim data
sim_data = errorbar(mt_lengths * site_size, avg_endtag_length, err_endtag_length, 'o', ...
    'MarkerSize', 4, 'LineWidth', 2);
sim_data.MarkerFaceColor = sim_data.MarkerEdgeColor;
% Plot exp data -- second plot is for horizontal error bars 
exp_data = errorbar(exp_mt_lengths, exp_endtag_lengths, ...
    exp_err_endtag_lengths, 'b^');
exp_data.MarkerFaceColor = exp_data.MarkerEdgeColor;
errorbarxy(exp_mt_lengths, exp_endtag_lengths, exp_err_mt_lengths, ...
    exp_err_endtag_lengths, {'b^', 'b', 'b'});



xlabel("Microtubule length (microns)", 'FontSize', 16);
ylabel("Endtag length (microns)", 'Fontsize', 16);
set(gca, 'FontSize', 16);
legend(["Experiment", "Simulation", "wut"], 'location', 'northwest'); 
%ylim([0 3]);
xlim([0 16]);
