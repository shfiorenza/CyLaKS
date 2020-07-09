base_name = "endtag";
mt_lengths = [250, 500, 750, 1000, 1250, 1750];
seeds = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
dir = "/home/shane/Projects/overlap_analysis/mgh_model";

n_mts = length(mt_lengths);
n_seeds = length(seeds);

endtag_lengths = zeros(n_mts, n_seeds);
for i_mt = 1 : n_mts
   for i_seed = 1 : n_seeds
      sim_name = sprintf("%s/%s_%i_%i", dir, base_name, mt_lengths(i_mt), seeds(i_seed));
      endtag_length = get_endtag_length(sim_name);
      endtag_lengths(i_mt, i_seed) = endtag_length;
   end
end

avg_endtags = zeros(n_mts, 1);
err_endtags = zeros(n_mts, 1);
for i_mt = 1 : n_mts
   avg_endtags(i_mt) = mean(endtag_lengths(i_mt, :));
   var = 0.0;
   for i_seed = 1 : n_seeds
       var = var + double(avg_endtags(i_mt) - endtag_lengths(i_mt, i_seed))^2 / (n_seeds - 1);
   end
   err_endtags(i_mt) = sqrt(var / n_seeds);
end

fig1 = figure();
set(fig1, 'Position', [50, 50, 2*720, 2*240])
errorbar(mt_lengths, avg_endtags, err_endtags);
