clear all;
applied_force = [0.05, 0.5, 5.0, 50.0]; % in pN
simNameBase = "mt_forceVel_" + applied_force;
seeds = [0, 1, 2, 3, 4, 5, 6, 7];
mt_lengths = [25, 125, 625, 3125, 15625]; % in n_sites
gammas = [4.83452e-06, 1.36262e-05, 4.74353e-05, 1.81917e-04, 7.37706e-04]; % in pN*s/site
site_size = 0.008; % in um
delta_t = 2.5e-5; % in s
n_steps = 4000000;
n_datapoints = 10000;
time_per_datapoint = delta_t * n_steps / n_datapoints;
n_sims = length(simNameBase);
n_seeds = length(seeds);
n_mts = length(mt_lengths);
expected_vel = zeros(n_sims, n_mts); % in um/s 
for i_sim = 1 : n_sims
    for i_mt = 1 : n_mts
        % Calculate velocity then convert to um/s
        expected_vel(i_sim, i_mt) = applied_force(i_sim) / gammas(i_mt) * site_size;
    end
end

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
       mt_raw_data = fread(mt_data_file, n_mts * n_datapoints, '*double');
       fclose(mt_data_file);
       mt_data = reshape(mt_raw_data, n_mts, n_datapoints);
       for i_mt = 1 : n_mts
           vel = gradient(double(mt_data(i_mt, :)) * site_size, time_per_datapoint);
           velocities(i_sim, i_seed, i_mt, :) = smoothdata(vel, 'movmean', 10);        
       end
   end
end

n_plot_points = 5;
start_point = 500;
data_increment = n_datapoints / n_plot_points;
end_point = start_point + (n_plot_points - 1) * data_increment;
avg_velocities = zeros(n_sims, n_mts, n_plot_points);
err_velocities = zeros(n_sims, n_mts, n_plot_points);
for i_sim = 1 : n_sims
    for i_mt = 1 : n_mts
        for i_entry = 1 : n_plot_points
            i_data = start_point + (i_entry - 1) * data_increment;
            avg_vel = mean(velocities(i_sim, :, i_mt, i_data));
            variance = 0.0;
            for i_seed = 1 : n_seeds
                diff_sq = (avg_vel - velocities(i_sim, i_seed, i_mt, i_data))^2;
                variance = variance + diff_sq/double(n_seeds - 1);
            end
            avg_velocities(i_sim, i_mt, i_entry) = avg_vel;
            err_velocities(i_sim, i_mt, i_entry) = sqrt(variance / n_seeds);
        end
    end
end

fig1 = figure();
set(fig1, 'Position', [50, 50, 1500, 500]);
t = linspace(start_point * time_per_datapoint, end_point * time_per_datapoint, n_plot_points);
for i_sim = 1 : n_sims
   subplot(1, n_sims + 1, i_sim)
   hold on
   for i_mt = 1 : n_mts
      errorbar(t, squeeze(avg_velocities(i_sim, i_mt, :)), squeeze(err_velocities(i_sim, i_mt, :)), ...
          '.', 'LineWidth', 2, 'MarkerSize', 14); 
   end
   for i_mt = 1 : n_mts
      plot([0 n_steps * delta_t], [expected_vel(i_sim, i_mt) expected_vel(i_sim, i_mt)], ...
          '--','Color', [0.25,0.25,0.25], 'LineWidth', 1.5); 
   end
   ax = gca;
   ax.FontSize = 12;   
   ylim([0 1.25*expected_vel(i_sim, 1)]);
   title({sprintf("F = %g pN", applied_force(i_sim)), " "}, 'FontSize', 12);
   xlabel("Time (s)", 'FontSize', 14);
   ylabel("Velocity (\mum/s)", 'FontSize', 14);

end
% Some trickery to make a legend common to all n_sims subplots
subplot(1, n_sims + 1, n_sims + 1);
hold on
for i_mt = 1 : n_mts
    errorbar([0,0],  [0,0], '.', 'LineWidth', 2, 'MarkerSize', 14);
end
plot(0,0, '--', 'Color', [0.25, 0.25, 0.25], 'LineWidth', 1.5);
ylim([1 2]);
xlim([1 2]);
axis off;
legendLabel = "Simulation (L = " + mt_lengths * site_size + " \mum)";
legendLabel(n_mts + 1) = "Theory (all lengths)";
legend(legendLabel, 'location', 'best', 'FontSize', 14);
sgtitle("Microtubule force-velocity relationship; thermal noise included", 'FontSize', 16);