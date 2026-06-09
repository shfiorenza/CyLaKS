sim_name = 'test_smol';

plot_histograms = false;
win_size = 600; % seconds
tip_pos = 0;
soma_pos = 40000;


% Load parameter structure
file_dir = '..';  % Default; only change if you move CyLaKS output files
params = load_parameters(sprintf('%s/%s', file_dir, sim_name));


filament_num_file = fopen(sprintf('%s/%s_filament_num.file', file_dir, sim_name));
filament_num_raw = fread(filament_num_file, '*double');
fclose(filament_num_file);
filament_lengths_file = fopen(sprintf('%s/%s_filament_lengths.file', file_dir, sim_name));
filament_lengths_raw = fread(filament_lengths_file, '*double');
fclose(filament_lengths_file);
filament_pos_file = fopen(sprintf('%s/%s_filament_pos.file', file_dir, sim_name));
filament_pos_raw = fread(filament_pos_file, '*double');
fclose(filament_pos_file);
filament_id_file = fopen(sprintf('%s/%s_filament_id.file', file_dir, sim_name));
filament_id_raw = fread(filament_id_file, '*double');
fclose(filament_id_file);


mt_num = filament_num_raw; 
mt_len = NaN(params.n_datapoints, max(mt_num));
mt_id = NaN(params.n_datapoints, max(mt_num));
mt_pos = NaN(params.n_datapoints, max(mt_num), 2, params.n_dims);

i_entry = 1;
j_entry = 1;
for i_datapoint = 1 : 1 : params.n_datapoints
    for i_mt = 1 : mt_num(i_datapoint)
        mt_len(i_datapoint, i_mt) = filament_lengths_raw(i_entry);
        mt_id(i_datapoint, i_mt) = filament_id_raw(i_entry);
        i_entry = i_entry + 1;
        for i_end = 1 : 1 : 2
            for i_dim = 1 : 1 : params.n_dims
                mt_pos(i_datapoint, i_mt, i_end, i_dim) = filament_pos_raw(j_entry);
                j_entry = j_entry + 1;
            end
        end
    end
end

n_mts_max = max(mt_num);
win_size_steps = win_size / params.time_per_datapoint;
vel_plus = nan(params.n_datapoints, n_mts_max);
vel_minus = nan(params.n_datapoints, n_mts_max);
n_plus_out = zeros(params.n_datapoints, 1);
n_minus_out = zeros(params.n_datapoints, 1);
for i_datapoint = win_size_steps + 1: 1 : params.n_datapoints
   for i_mt = 1 : mt_num(i_datapoint)
        id = mt_id(i_datapoint, i_mt);
        cur_pos_plus = mt_pos(i_datapoint, i_mt, 1, 1);
        cur_pos_minus = mt_pos(i_datapoint, i_mt, 2, 1);
        i_old = find(mt_id(i_datapoint - win_size_steps, :) == id, 1);
        if ~isempty(i_old) %&& cur_pos_plus < 10000
            old_pos_plus = mt_pos(i_datapoint - win_size_steps, i_old, 1, 1);
            old_pos_minus = mt_pos(i_datapoint - win_size_steps, i_old, 2, 1);
            vp = (cur_pos_plus - old_pos_plus) / win_size  * 60 / 1000;
            vm = (cur_pos_minus - old_pos_minus) / win_size * 60 / 1000;
            vel_plus(i_datapoint, i_mt) = vp;
            vel_minus(i_datapoint, i_mt) = vm;
            if cur_pos_plus > cur_pos_minus
                n_minus_out(i_datapoint) = n_minus_out(i_datapoint) + 1;
            else
                n_plus_out(i_datapoint) = n_plus_out(i_datapoint) + 1;
            end
        end
   end
end

if plot_histograms
    fig = figure('Position', [50 50 720 360]);
end

avg_vel_plus = NaN(params.n_datapoints, 1);
avg_vel_minus = NaN(params.n_datapoints, 1);
err_vel_minus = NaN(params.n_datapoints, 1);

for i_datapoint = win_size_steps + 1: 1 : params.n_datapoints
    serial_plus = rmmissing(reshape(vel_plus(i_datapoint, :), [], 1));
    avg_plus = mean(serial_plus);
    sd_plus = std(serial_plus);
    sem_plus = sd_plus / sqrt(length(serial_plus));

    serial_minus = rmmissing(reshape(vel_minus(i_datapoint, :), [], 1));
    avg_minus = mean(serial_minus); 
    sd_minus = std(serial_minus);
    sem_minus = sd_minus / sqrt(length(serial_minus));

    avg_vel_plus(i_datapoint) = avg_plus;
    avg_vel_minus(i_datapoint) = avg_minus;
    err_vel_minus(i_datapoint) = sem_minus;

    if ~plot_histograms
        continue;
    end
    
    clf;  
    histogram(squeeze(vel_plus(i_datapoint, :)), 'FaceColor',[1 0 0], 'FaceAlpha', 0.5, 'EdgeColor', 'none')
    hold on
    histogram(squeeze(vel_minus(i_datapoint, :)), 'FaceColor',[0 0 1], 'FaceAlpha', 0.5, 'EdgeColor', 'none')
    axis square
    dim = [0.75 0.6 .3 .3];
    time = (i_datapoint) * params.time_per_datapoint;
    str = sprintf('Time: %#.2f seconds', time);
    annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on');
    dim3 = [0.75 0.55 0.2 0.2];
    str3 = sprintf('Plus: %#.2f +/- %#.2f um/min', avg_plus, sem_plus);
    annotation('textbox', dim3, 'String', str3, 'FitBoxToText', 'on');
    dim4 = [0.75 0.4 0.2 0.2];
    str4 = sprintf('Minus: %#.2f +/- %#.2f um/min', avg_minus, sem_minus);
    annotation('textbox', dim4, 'String', str4, 'FitBoxToText', 'on');
drawnow()
end

polarity = n_plus_out ./ (n_plus_out + n_minus_out);

fig = figure('Position', [0 50 720 360]);
%plot(linspace(0, params.t_run/60, params.n_datapoints), smooth(avg_vel_plus), 'LineWidth', 3);
hold on
yyaxis left
plot(linspace(0, params.t_run/60, params.n_datapoints), avg_vel_minus, 'LineWidth', 3);
ylabel("Velocity (um/min)");
yyaxis right
plot(linspace(0, params.t_run/60, params.n_datapoints), polarity, 'LineWidth', 3);
ylabel("Polarity (unitless)");
%errorbar(linspace(0, params.t_run/60, params.n_datapoints), avg_vel_minus, err_vel_minus, 'LineWidth', 3);
xlim([0 params.t_run / 60])
xlabel("Time (min)");

return

active_ids = zeros(n_mts_max);
last_pos = zeros(n_mts_max);

%pos = zeros(params.n_datapoints, 10);
minus_pos = NaN(params.n_datapoints, max(mt_num));
for i_datapoint = 2 : 1 : params.n_datapoints
   for i_mt = 1 : mt_num(i_datapoint)
    active_ids(i_mt) = mt_id(i_datapoint, i_mt);
    last_pos(i_mt) = minus_pos(i_datapoint - 1, i_mt, 2, 1);


      
     minus_pos(i_datapoint, i_mt) = mt_pos(i_datapoint, i_mt, 2, 1);
     pos(i_datapoint, i_mt) = mt_pos(i_datapoint, i_mt, 2, 1);
   end
end

fig = figure();
%plot(minus_pos);
%hold on
avg = mean(minus_pos, 2, 'omitnan');
plot(linspace(0, params.t_run / 60, params.n_datapoints), avg / 1000, 'LineWidth', 3);
fontname("arial");
fontsize(14, "points");
xlabel("Time (min)", "FontSize", 18);
ylabel("Avg. minus-end position (um)", "FontSize", 18);
ylim([0 40]);
grid on
box off


fig = figure();

plot(linspace(0, params.t_run / 60, params.n_datapoints), mt_num, 'LineWidth', 3)
fontname("arial");
fontsize(14, "points");
xlabel("Time (min)", "FontSize", 18);
ylabel("Microtubule count", "FontSize", 18);
ylim([0 700]);
