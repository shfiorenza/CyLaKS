%{
clear variables;
file_dir = '/home/shane/projects/CyLaKS/%s';
sim_name_base = ["run_mt_forceVel/mt_forceVel"];
%sim_name_base = ["mt_forceVel0.5", "mt_forceVel_5.0", "mt_forceVel_50.0"];
seeds = [0, 1, 2, 3, 4, 5]; 
applied_force = [1 10:10:100];
%applied_force = [1, 5.0, 10, 50, 100]; % in pN
labels = ["1"];
for i = 10 : 10 : 40
   labels = [labels num2str(i)]; 
end
labels = [labels "50.0"];
for i = 60 : 10 : 100
   labels = [labels num2str(i)]; 
end

n_dims = 2;
% Open log file and parse it into param labels & their values
log_file = sprintf(file_dir, sprintf('%s.log', sim_name_base));
if(~isempty(seeds) && ~isempty(applied_force))
    log_file = sprintf(file_dir, sprintf('%s_%s_%i.log', sim_name_base, labels(1), seeds(1)));
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
n_sims = length(applied_force); %length(sim_name_base);
expected_vel_par = zeros(n_sims, n_mts); % in um/s
expected_vel_perp = zeros(n_sims, n_mts); % in um/s
for i_sim = 1:n_sims
    for i_mt = 1:n_mts
        % Calculate velocity then convert to um/s
        expected_vel_par(i_sim, i_mt) = (applied_force(i_sim) / gamma_par(i_mt));
        expected_vel_perp(i_sim, i_mt) = (applied_force(i_sim) / gamma_perp(i_mt));
    end
end

velocities_par = zeros(n_sims, n_seeds, n_mts, n_datapoints);
velocities_perp = zeros(n_sims, n_seeds, n_mts, n_datapoints);
for i_sim = 1 : n_sims
    for i_seed = 1 : n_seeds
        sim_name = sim_name_base;
        if(~isempty(seeds))
            sim_name = sprintf('%s_%s_%i', sim_name_base,labels(i_sim), seeds(i_seed));
        end
        filename = sprintf(file_dir, sprintf('%s_filament_pos.file', sim_name))
        file = fopen(filename);
        data = fread(file, 2*n_dims * n_mts * n_datapoints, '*double');
        filament_pos = reshape(data, n_dims, 2, n_mts, n_datapoints);
        fclose(file);

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
            velocities_par(i_sim, i_seed, i_mt, :) = vel_x; %smoothdata(vel_x, 'movmean', 10);
            velocities_perp(i_sim, i_seed, i_mt, :) = vel_y;
        end
    end
end
%}
%scan_window = 10;
%i_win = scan_window / time_per_datapoint;
avg_velocities_par = zeros(n_sims, n_mts);
err_velocities_par = zeros(n_sims, n_mts);
for i_sim = 1 : n_sims
    for i_mt = 1 : n_mts
        avg_vel = mean(velocities_par(i_sim, :, i_mt, :), [2 4]);
        variance = 0.0;
        for i_seed = 1 : n_seeds
            diff_sq = (avg_vel - mean(velocities_par(i_sim, i_seed, i_mt, :)))^2;
            variance = variance + diff_sq / double(n_seeds - 1);
        end
        avg_velocities_par(i_sim, i_mt) = avg_vel;
        err_velocities_par(i_sim, i_mt) = sqrt(variance);
    end
end
avg_velocities_perp = zeros(n_sims, n_mts);
err_velocities_perp = zeros(n_sims, n_mts);
for i_sim = 1 : n_sims
    for i_mt = 1 : n_mts
        avg_vel = mean(velocities_perp(i_sim, :, i_mt, :), [2 4]);
        variance = 0.0;
        for i_seed = 1 : n_seeds
            diff_sq = (avg_vel - mean(velocities_perp(i_sim, i_seed, i_mt, :)))^2;
            variance = variance + diff_sq / double(n_seeds - 1);
        end
        avg_velocities_perp(i_sim, i_mt) = avg_vel;
        err_velocities_perp(i_sim, i_mt) = sqrt(variance);
    end
end
%}

data_color = [0, 0.4470, 0.7410; 0.9290, 0.6940, 0.1250; 0.4660, 0.6740, 0.1880];
line_color = [0.8500, 0.3250, 0.0980;0.4940, 0.1840, 0.5560; 0.3010, 0.7450, 0.9330];
color = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980;  0.9290, 0.6940, 0.1250];

fig1 = figure();
set(fig1, 'Position', [50, 50, 720, 540]);
hold all;

for i_mt = 1 : n_mts
    %sim_data = errorbar(applied_force, squeeze(0.001 * avg_velocities_par(:, i_mt)), squeeze(0.001 * err_velocities_par(:, i_mt)), ...
    sim_data = plot(applied_force, squeeze(0.001 * avg_velocities_par(:, i_mt)), ...
        'o', 'MarkerSize', 12, 'MarkerEdgeColor', color(i_mt, :));
    sim_data.MarkerFaceColor = sim_data.MarkerEdgeColor;
    sim_data.Color = sim_data.MarkerEdgeColor;
    %sim_data = errorbar(applied_force, squeeze(0.001 * avg_velocities_perp(:, i_mt)), squeeze(0.001 * err_velocities_perp(:, i_mt)), ...
    sim_data = plot(applied_force, squeeze(0.001 * avg_velocities_perp(:, i_mt)), ...
        'sq', 'MarkerSize', 12, 'MarkerEdgeColor', color(i_mt, :));
    sim_data.MarkerFaceColor = sim_data.MarkerEdgeColor;
    sim_data.Color = sim_data.MarkerEdgeColor;
    plot([0 100], [0 0.1/gamma_par(i_mt)], 'LineWidth', 2, 'Color', [0.6 0.6 0.6]);
    plot([0 100], [0 0.1/gamma_perp(i_mt)], 'LineWidth', 2, 'Color', [0.6 0.6 0.6]);
end
xlim([-3 105]);
ylim([-30 750]);
xticks([0 25 50 75 100]);
yticks([0 150 300 450 600]);
set(gca, 'FontSize', 22);
xlabel("Applied force (pN)", 'FontSize', 22);
ylabel("Velocity (\mum/s)", 'FontSize', 22);
legendLabel = "L = " + num2str(n_sites' * 0.008) + " \mum";
legendLabel(n_mts + 1) = "Theory (all lengths)";
h = get(gca,'Children');
n_entries = length(h);
h_array = [];
for i = n_entries : - 1 : 1
      h_array = [h_array h(i)];
end
set(gca,'Children',h_array)
legend(h([n_entries, n_entries-4, n_entries-8, n_entries-2, n_entries-1, n_entries-5, n_entries-9]), ... 
    ["L = 1 um (par)", "L = 5 um (par)", "L = 20 um (par)", "Theory", "L = 1um (perp)", "L = 5 um (perp)", "L = 20 um (perp)"], ...
    'location', 'northwest', 'FontSize', 18, 'NumColumns', 2);
legend('boxoff') 
%{
h = get(gca,'Children');
legend([h(5) h(3) h(1)],legendLabel, 'location', 'northwest', 'FontSize', 20);
legend('boxoff');
%}
%set(gca, 'XScale', 'log');
