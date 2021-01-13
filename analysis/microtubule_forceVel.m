clear variables;
file_dir = '/home/shane/projects/CyLaKS/%s';
sim_name_base = ['test'];
%sim_name_base = ["mt_forceVel_0.5", "mt_forceVel_5.0", "mt_forceVel_50.0"];
seeds = []; %[0, 1, 2, 3, 4]; 
applied_force = [0.5, 5.0, 50.0]; % in pN

n_dims = 2;
% Open log file and parse it into param labels & their values
log_file = sprintf(file_dir, sprintf('%s.log', sim_name_base));
if(~isempty(seeds))
    log_file = sprintf(file_dir, sprintf('%s_%i.log', sim_name_base, seeds(0)));
end
log = textscan(fileread(log_file), '%s %s', 'Delimiter', '=');
params = log{1, 1};
values = log{1, 2};
% Read in number of MTs
n_mts = sscanf(values{contains(params, 'count ')}, '%g');
n_sites = zeros(1, n_mts);
ell = zeros(1, n_mts);
gamma_par = zeros(1, n_mts);
gamma_perp = zeros(1, n_mts);
for i_mt = 1 : n_mts
    string = sprintf('n_sites[%i] ', i_mt - 1);
    n_sites(i_mt) = sscanf(values{contains(params, string)}, '%i');
    string = sprintf('length[%i] ', i_mt - 1);
    ell(i_mt) = sscanf(values{contains(params, string)}, '%g');
    string = sprintf('gamma_par[%i] ', i_mt - 1);
    gamma_par(i_mt) = sscanf(values{contains(params, string)}, '%g');
    string = sprintf('gamma_perp[%i] ', i_mt - 1);
    gamma_perp(i_mt) = sscanf(values{contains(params, string)}, '%g');
end
% Read in system params
dt = sscanf(values{contains(params, 'dt ')}, '%g');
steps_per_datapoint = str2double(values{contains(params, 'n_steps_per_snapshot ')});
time_per_datapoint = dt * steps_per_datapoint;
n_datapoints = str2double(values{contains(params, 'n_datapoints ')});
% Use actual recorded number of datapoints to parse thru data/etc
if any(contains(params, 'N_DATAPOINTS ') ~= 0)
    n_datapoints = str2double(values{contains(params, 'N_DATAPOINTS ')});
end
n_seeds = 1;
if(~isempty(seeds))
    n_seeds = length(seeds);
end
n_sims = length(sim_name_base);
expected_vel = zeros(n_sims, n_mts); % in um/s
for i_sim = 1:n_sims
    for i_mt = 1:n_mts
        % Calculate velocity then convert to um/s
        expected_vel(i_sim, i_mt) = (applied_force(i_sim) / gammas(i_mt)) / 1000;
    end
end

velocities = zeros(n_sims, n_seeds, n_mts, n_datapoints);
for i_sim = 1 : n_sims
    for i_seed = 1 : n_seeds
        sim_name = sim_name_base;
        if(~isempty(seeds))
            sim_name = sprintf('%s_%i', sim_name_base, seeds(i_seed));
        end
        filename = '%s_filament_pos.file';
        file = fopen(sprintf(file_dir, sprintf(filename, sim_name)));
        data = fread(file, 2*n_dims * n_mts * n_datapoints, '*double');
        filament_pos = reshape(data, n_dims, 2, n_mts, n_datapoints);
        fclose(file);

        mt_data = reshape(mt_raw_data, n_mts, n_datapoints);

        for i_mt = 1:n_mts
            com_x = zeros(1, n_datapoints);
            com_y = zeros(1, n_datapoints); 
            for i_data = 1 : n_datapoints
                plus_pos = filament_pos(:, 1, i_mt, i_data);
                minus_pos = filament_pos(:, 2, i_mt, i_data);
                com_x(i_data) = double(plus_pos(1) + minus_pos(1))/2.0;
                com_y(i_data) = double(plus_pos(2) + minus_pos(2))/2.0;
            end
            vel_x = gradient(com_x, time_per_datapoint);
            vel_y = gradient(com_y, time_per_datapoint);
            velocities(i_sim, i_seed, i_mt, :) = smoothdata(vel_x, 'movmean', 10);
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
                variance = variance + diff_sq / double(n_seeds - 1);
            end
            avg_velocities(i_sim, i_mt, i_entry) = avg_vel;
            err_velocities(i_sim, i_mt, i_entry) = sqrt(variance);
        end
    end
end

fig1 = figure();
set(fig1, 'Position', [50, 50, 1500, 500]);
t = linspace(0, n_steps * delta_t, n_plot_points);
for i_sim = 1:n_sims
    subplot(1, n_sims + 1, i_sim)
    hold on

    for i_mt = 1:n_mts
        errorbar(t, squeeze(avg_velocities(i_sim, i_mt, :)), squeeze(err_velocities(i_sim, i_mt, :)), ...
            'o', 'LineWidth', 2, 'MarkerSize', 10);
    end

    for i_mt = 1:n_mts
        plot([0 n_steps * delta_t], [expected_vel(i_sim, i_mt) expected_vel(i_sim, i_mt)], ...
            '--', 'Color', [0.25, 0.25, 0.25], 'LineWidth', 1.5);
    end

    ax = gca;
    ax.FontSize = 12;
    ylim([0 1.1 * expected_vel(i_sim, 1)]);
    title({sprintf("Force = %g pN", applied_force(i_sim)), " "}, 'FontSize', 12);
    xlabel("Time (s)", 'FontSize', 12);
    ylabel("Velocity (\mum/s)", 'FontSize', 12);

end

% Some trickery to make a legend common to all n_sims subplots
subplot(1, n_sims + 1, n_sims + 1);
hold on

for i_mt = 1:n_mts
    errorbar([0, 0], [0, 0], 'o', 'LineWidth', 2, 'MarkerSize', 10);
end

plot(0, 0, '--', 'Color', [0.25, 0.25, 0.25], 'LineWidth', 1.5);
ylim([1 2]);
xlim([1 2]);
axis off;
legendLabel = "Sim data (length = " + num2str(mt_lengths' * site_size) + " \mum)";
legendLabel(n_mts + 1) = "Theory (all lengths)";
legend(legendLabel, 'location', 'northwestoutside', 'FontSize', 10);
