clear all;
simNameBase = ["mt_forceVel_0.5", "mt_forceVel_5.0", "mt_forceVel_50.0"];
applied_force = [0.5, 5.0, 50.0]; % in pN
seeds = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
mt_lengths = [25, 125, 625, 3125]; % in n_sites
%gammas = [3.02e-06, 1.21e-05, 4.83e-05]; % in (pN*s)/nm
gammas = [6.04315e-07, 1.70328e-06, 5.92942e-06, 2.27396e-05];
n_sims = length(simNameBase);
n_seeds = length(seeds);
n_mts = length(mt_lengths);
expected_vel = zeros(n_sims, n_mts); % in um/s 
for i_sim = 1 : n_sims
    for i_mt = 1 : n_mts
        % Calculate velocity then convert to um/s
        expected_vel(i_sim, i_mt) = (applied_force(i_sim) / gammas(i_mt)) / 1000;
    end
end
site_size = 0.008; % in um
delta_t = 2.5e-5; % in s
n_steps = 6000000;
n_datapoints = 10000;
time_per_datapoint = delta_t * n_steps / n_datapoints;

fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
mtFileStruct = '%s_mt_coord.file';

velocities = zeros(n_sims, n_seeds, n_mts, n_datapoints);

for i_sim = 1 : n_sims
   for i_seed = 1 : n_seeds
       simName = simNameBase(i_sim);
       if(n_seeds > 1)
          simName = sprintf("%s_%i", simNameBase(i_sim), seeds(i_seed)); 
       end
       mtFileName = sprintf(fileDirectory, sprintf(mtFileStruct, simName));
       mt_data_file = fopen(mtFileName);
       mt_raw_data = fread(mt_data_file, n_mts * n_datapoints, '*int');
       fclose(mt_data_file);
       mt_data = reshape(mt_raw_data, n_mts, n_datapoints);
       for i_mt = 1 : n_mts
           vel = gradient(double(mt_data(i_mt, :)) * site_size, time_per_datapoint);
           velocities(i_sim, i_seed, i_mt, :) = smoothdata(vel, 'movmean', 10);        
       end
   end
end

n_plot_points = 10;
data_increment = n_datapoints / n_plot_points;
avg_velocities = zeros(n_sims, n_mts, n_plot_points);
err_velocities = zeros(n_sims, n_mts, n_plot_points);
for i_sim = 1 : n_sims
    for i_mt = 1 : n_mts
        for i_entry = 1 : n_plot_points
            i_data = 1 + (i_entry - 1) * data_increment;
            avg_vel = mean(velocities(i_sim, :, i_mt, i_data));
            variance = 0.0;
            for i_seed = 1 : n_seeds
                diff_sq = (avg_vel - velocities(i_sim, i_seed, i_mt, i_data))^2;
                variance = variance + diff_sq/double(n_seeds - 1);
            end
            avg_velocities(i_sim, i_mt, i_entry) = avg_vel;
            err_velocities(i_sim, i_mt, i_entry) = sqrt(variance);
        end
    end
end

fig1 = figure();
set(fig1, 'Position', [50, 50, 1500, 500]);
t = linspace(0, n_steps * delta_t, n_plot_points);
for i_sim = 1 : n_sims
   subplot(1, n_sims + 1, i_sim)
   hold on
   for i_mt = 1 : n_mts
      errorbar(t, squeeze(avg_velocities(i_sim, i_mt, :)), squeeze(err_velocities(i_sim, i_mt, :)), ...
          'o', 'LineWidth', 2, 'MarkerSize', 10); 
   end
   for i_mt = 1 : n_mts
      plot([0 n_steps * delta_t], [expected_vel(i_sim, i_mt) expected_vel(i_sim, i_mt)], ...
          '--','Color', [0.25,0.25,0.25], 'LineWidth', 1.5); 
   end
   ax = gca;
   ax.FontSize = 12;   
   ylim([0 1.1*expected_vel(i_sim, 1)]);
   title({sprintf("Force = %g pN", applied_force(i_sim)), " "}, 'FontSize', 12);
   xlabel("Time (s)", 'FontSize', 12);
   ylabel("Velocity (\mum/s)", 'FontSize', 12);

end
% Some trickery to make a legend common to all n_sims subplots
subplot(1, n_sims + 1, n_sims + 1);
hold on
for i_mt = 1 : n_mts
errorbar([0,0],  [0,0], 'o', 'LineWidth', 2, 'MarkerSize', 10);
end
plot(0,0, '--', 'Color', [0.25, 0.25, 0.25], 'LineWidth', 1.5);
ylim([1 2]);
xlim([1 2]);
axis off;
legendLabel = "Sim data (length = " + num2str(mt_lengths' * site_size) + " \mum)";
legendLabel(n_mts + 1) = "Theory (all lengths)";
legend(legendLabel, 'location', 'northwestoutside', 'FontSize', 10);

