clear variables;

sim_name = 'test_long_noFric_noGrowth_4xLength';
sim_name = 'test_longID_new_0.01xFric_15xF0';

%output_movie_name = 'mov_all_test_longID_new_0.01xFric_1xF0'; %'mov_all_test_long_noFric_noGrowth_4xLength';
output_movie_name = 'mov_all_posMTRF_longID_new_0.01xFric_15xF0'; %'mov_all_test_long_noFric_noGrowth_4xLength';

win_size = 600; % seconds
%axon_tip_pos = 10000;
soma_pos = 0;

site_size = 8.2; %nm

size = [50 50 1200 500];

start_point = 1;
end_point = -1;  % set to -1 to run until end of data

datapoints_per_plot = 10;
movie_duration = 15; % in seconds

% Load parameter structure
file_dir = '..';  % Default; only change if you move CyLaKS output files
params = load_parameters(sprintf('%s/%s', file_dir, sim_name));

%params.n_datapoints = 1320; %1000; %24;
%params.t_run = 13200; %10000;

if end_point == -1
    end_point = params.n_datapoints;
end
n_plots = (end_point - start_point + 1) / datapoints_per_plot;
r_prot = (params.site_size*1000);

% Initialize videowriter object
v = VideoWriter(output_movie_name, 'MPEG-4');
v.FrameRate = n_plots / movie_duration;
open(v);
frame_box = size;

filament_num_file = fopen(sprintf('%s/%s_filament_num.file', file_dir, sim_name));
filament_num_raw = fread(filament_num_file, '*double');
fclose(filament_num_file);
filament_lengths_file = fopen(sprintf('%s/%s_filament_lengths.file', file_dir, sim_name));
filament_lengths_raw = fread(filament_lengths_file, '*double');
fclose(filament_lengths_file);
filament_pos_file = fopen(sprintf('%s/%s_filament_pos.file', file_dir, sim_name));
filament_pos_raw = fread(filament_pos_file, '*double');
fclose(filament_pos_file);
filament_tip_file = fopen(sprintf('%s/%s_filament_pos_tip.file', file_dir, sim_name));
filament_tip_raw = fread(filament_tip_file, '*double');
fclose(filament_tip_file);
filament_id_file = fopen(sprintf('%s/%s_filament_id.file', file_dir, sim_name));
filament_id_raw = fread(filament_id_file, '*double');
fclose(filament_id_file);

axon_tip_pos_file = fopen(sprintf('%s/%s_axon_tip_pos.file', file_dir, sim_name));
axon_tip_pos_raw = fread(axon_tip_pos_file, '*double');
fclose(axon_tip_pos_file);
axon_tip_force_file = fopen(sprintf('%s/%s_axon_tip_force.file', file_dir, sim_name));
axon_tip_force_raw = fread(axon_tip_force_file, '*double');
fclose(axon_tip_force_file);
%}

mt_num = filament_num_raw;
axon_tip_pos = axon_tip_pos_raw; %40000 * ones(params.n_datapoints, 1); %axon_tip_pos_raw;
axon_tip_force = axon_tip_force_raw;
mt_len = NaN(params.n_datapoints, max(mt_num));
mt_pos = NaN(params.n_datapoints, max(mt_num), 2, params.n_dims);
mt_id = NaN(params.n_datapoints, max(mt_num));
mt_tip = NaN(params.n_datapoints, max(mt_num), params.n_dims);

i_entry = 1;
j_entry = 1;
k_entry = 1;
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
        for i_dim = 1 : 1 : params.n_dims
            mt_tip(i_datapoint, i_mt, i_dim) = filament_tip_raw(k_entry);
            k_entry = k_entry + 1;
        end
    end
end

% Open figure and set to desired size (each frame must be this same size)
fig1 = figure('Position', size);

n_mts_max = max(mt_num);
win_size_steps = win_size / params.time_per_datapoint;
vel_plus = nan(params.n_datapoints, n_mts_max, 3);
vel_plusStable = nan(params.n_datapoints, n_mts_max, 3);
vel_minus = nan(params.n_datapoints, n_mts_max, 3);
n_plus_out = zeros(params.n_datapoints, 3);
n_minus_out = zeros(params.n_datapoints, 3);
for i_datapoint = win_size_steps + 1: 1 : params.n_datapoints
    for i_mt = 1 : mt_num(i_datapoint)
        id = mt_id(i_datapoint, i_mt);
        cur_pos_plus = mt_pos(i_datapoint, i_mt, 1, 1);
        cur_pos_plusStable = mt_tip(i_datapoint, i_mt, 1);
        cur_pos_minus = mt_pos(i_datapoint, i_mt, 2, 1);
        i_old = find(mt_id(i_datapoint - win_size_steps, :) == id, 1);
        if ~isempty(i_old)
            old_pos_plus = mt_pos(i_datapoint - win_size_steps, i_old, 1, 1);
            old_pos_plusStable = mt_tip(i_datapoint - win_size_steps, i_old, 1);
            old_pos_minus = mt_pos(i_datapoint - win_size_steps, i_old, 2, 1);
            vp = (cur_pos_plus - old_pos_plus) / win_size  * 60 / 1000;
            vps = (cur_pos_plusStable - old_pos_plusStable) / win_size * 60 / 1000;
            vm = (cur_pos_minus - old_pos_minus) / win_size * 60 / 1000;
            if cur_pos_plusStable > (axon_tip_pos(i_datapoint) - soma_pos) * 2 / 3.0
                i_region = 1;
            elseif cur_pos_plusStable > (axon_tip_pos(i_datapoint) - soma_pos) / 3.0
                i_region = 2;
            else
                i_region = 3;
            end
            vel_plus(i_datapoint, i_mt, i_region) = vp;
            vel_plusStable(i_datapoint, i_mt, i_region) = vps;
            vel_minus(i_datapoint, i_mt, i_region) = vm;
            if cur_pos_plus < cur_pos_minus
                n_minus_out(i_datapoint, i_region) = n_minus_out(i_datapoint, i_region) + 1;
            else
                n_plus_out(i_datapoint, i_region) = n_plus_out(i_datapoint, i_region) + 1;
            end
        end
    end
end


%min_x_global = min(mt_pos(:, :, :, 1), [], "all");
%max_x_global = max(mt_pos(:, :, :, 1), [], "all");
%span = max_x_global - min_x_global;
%n_sites = span / 8.2;
%occu_data_stable = zeros(n_plots, ceil(n_sites));
%occu_data_unstable = zeros(n_plots, ceil(n_sites));
avg_vel_minus = NaN(n_plots, 3);
avg_polarity = NaN(n_plots, 3);
% Run through all datapoints; each one is a frame in our movie
i_plot = 0;
for i_datapoint = start_point : datapoints_per_plot : end_point
    i_plot = i_plot + 1;
    min_x = min(mt_pos(i_datapoint:i_datapoint+datapoints_per_plot-1, :, :, 1), [], "all");
    max_x = max(mt_pos(i_datapoint:i_datapoint+datapoints_per_plot-1, :, :, 1), [], "all");
    if min_x ~= min_x || max_x ~= max_x
        return
    end
    % Clear figure so that it only displays figures from current datapoint
    clf;

    % Set Axes properties
    ax = subplot(3, 3, 1:3);
    
    %ax = axes('Units', 'normalized', 'Position', [0.075 0.15 0.9 0.8]);
    %set(gca,'xdir','reverse');%,'ydir','reverse')
    hold all;
    %ax.YLim = [-40 100];
    xlabel("Position (microns)");
    ylabel("Position (nanometers)");
    %xlim([min_x_global - 1000 max_x_global + 1000])
    %xlim([min_x - 250 max_x + 250])
    %xlim([-250 max_x + 250])
    xlim([(soma_pos-250)/1000.0 (axon_tip_pos(i_datapoint) + 250)/1000.0])
    ylim([-15 35]);
    %xticks([]);
    %ax.YLabel.String = 'y position (nm)';
    for i_mt = 1:1:mt_num(i_datapoint)
        plus_pos = mt_pos(i_datapoint, i_mt, 1, :);
        minus_pos = mt_pos(i_datapoint, i_mt, 2, :);
        tip_pos = mt_tip(i_datapoint, i_mt, :);
        if plus_pos(1) > minus_pos(1)
            polarity = 0;
            color = [0.25 0.25 0.25];
        else
            polarity = 1;
            color = [0.7 0.7 0.7];
        end
        line([plus_pos(1)-r_prot/2, minus_pos(1)-r_prot/2]/1000.0,[plus_pos(2), minus_pos(2)], ...
            'LineWidth', 2, 'Color', color);
        line([plus_pos(1)-r_prot/2, tip_pos(1)-r_prot/2]/1000.0,[plus_pos(2), tip_pos(2)], ...
            'LineWidth', 2, 'Color', [0 1 0]);
    end

    %}
    %ax = axes('Units', 'normalized', 'Position', [0.075 0.17 0.9 0.8]);
    ax = subplot(3, 3, 4:6);

    span = max_x - min_x;
    n_sites = span / site_size;
    occu_data_stable = zeros(1, ceil(n_sites));
    occu_data_unstable = zeros(1, ceil(n_sites));
    %set(gca,'xdir','reverse');%,'ydir','reverse')
    hold all;
    for i_entry = i_datapoint : i_datapoint + datapoints_per_plot - 1
        for i_mt = 1 : mt_num(i_entry)
            plus_pos = mt_pos(i_entry, i_mt, 1, :);
            minus_pos = mt_pos(i_entry, i_mt, 2, :);
            tip_pos = mt_tip(i_entry, i_mt, :);
            len_sites = mt_len(i_entry, i_mt);
            len_tip = abs(tip_pos(1) - plus_pos(1));
            pos_start = min(mt_pos(i_entry, i_mt, :, 1));
            i_start = max(0, floor((pos_start - min_x)/site_size)) + 1;
            %disp(len_tip)
            
            for i = i_start : i_start + len_sites - 2
                if len_tip == 0
                    occu_data_stable(1, i) = occu_data_stable(1, i) + 1/datapoints_per_plot;
                else
                    occu_data_unstable(1, i) = occu_data_unstable(1, i) + 1/datapoints_per_plot;
                end
            end
            %}
        end
    end
    plot(linspace(min_x/1000.0, max_x/1000.0, ceil(n_sites)), occu_data_stable(1, :), 'LineWidth', 3, 'Color', 'k')
    hold on
    plot(linspace(min_x/1000.0, max_x/1000.0, ceil(n_sites)), occu_data_unstable(1, :), 'LineWidth', 3, 'Color', [0 1 0])
    legendLabel = {"Stable MTs", "Unstable MTs"};
    legend(legendLabel, 'location', 'northeast');
    ylabel("Microtubule density (A.U.)")
    %xlim([min_x_global/1000.0 - 1 max_x_global/1000.0 + 1])
    %xlim([(min_x - 250)/1000.00 (max_x + 250)/1000.0])
    %xlim([(-250)/1000.00 (max_x + 250)/1000.0])
    xlim([(soma_pos-250)/1000.00 (axon_tip_pos(i_datapoint) + 250)/1000.0])



    subplot(3, 3, 7); %nexttile([1 3])
    title ("Near soma");
    hold on
    yyaxis left
    serial_minus = rmmissing(reshape(vel_minus(i_datapoint, :, 3), [], 1));
    %serial_minus = rmmissing(reshape(vel_plusStable(i_frame, :, 3), [], 1));
    avg_minus = mean(serial_minus);
    avg_vel_minus(i_plot, 3) = avg_minus;
    avg_polarity(i_plot, 3) = n_plus_out(i_datapoint, 3) / (n_plus_out(i_datapoint, 3) + n_minus_out(i_datapoint, 3));
    plot(linspace(0, params.t_run/60, n_plots), -1*avg_vel_minus(:, 3), 'LineWidth', 2)
    xlim([0 params.t_run/60]);
    yline(0, '--');
    ylabel("MT-RF (um/min)");
    ylim([-0.5 0.5]);
    %ax.YAxis(1).Color = [128 0 128]/255;
    yyaxis right
    plot(linspace(0, params.t_run/60, n_plots), avg_polarity(:, 3), 'LineWidth', 2)
    ylim([0 1.1]);

    subplot(3, 3, 8); %nexttile([1 3])
    title("Midzone");
    hold on
    yyaxis left
    serial_minus = rmmissing(reshape(vel_minus(i_datapoint, :, 2), [], 1));
    %serial_minus = rmmissing(reshape(vel_plusStable(i_frame, :, 2), [], 1));
    avg_minus = mean(serial_minus);
    avg_vel_minus(i_plot, 2) = avg_minus;
    avg_polarity(i_plot, 2) = n_plus_out(i_datapoint, 2) / (n_plus_out(i_datapoint, 2) + n_minus_out(i_datapoint, 2));
    plot(linspace(0, params.t_run/60, n_plots), -1*avg_vel_minus(:, 2), 'LineWidth', 2)
    xlim([0 params.t_run/60]);
    %ax.YAxis(1).Color = [128 0 128]/255;
    yline(0, '--');
    ylim([-0.5 0.5]);
    yyaxis right
    plot(linspace(0, params.t_run/60, n_plots), avg_polarity(:, 2), 'LineWidth', 2)
    ylim([0 1.1]);
    xlabel("Time (minutes)");

        subplot(3, 3, 9); %nexttile([1 3])
        title("Near tip");
    hold on
    yyaxis left
    serial_minus = rmmissing(reshape(vel_minus(i_datapoint, :, 1), [], 1));
    %serial_minus = rmmissing(reshape(vel_plusStable(i_frame, :, 1), [], 1));
    avg_minus = mean(serial_minus);
    avg_vel_minus(i_plot, 1) = avg_minus;
    avg_polarity(i_plot, 1) = n_plus_out(i_datapoint, 1) / (n_plus_out(i_datapoint, 1) + n_minus_out(i_datapoint, 1));
    plot(linspace(0, params.t_run/60, n_plots), -1*avg_vel_minus(:, 1), 'LineWidth', 2)
    yline(0, '--');
    %ax.YAxis(1).Color = [128 0 128]/255;
    xlim([0 params.t_run/60]);
    ylim([-0.5 0.5]);
    yyaxis right
    plot(linspace(0, params.t_run/60, n_plots), avg_polarity(:, 1), 'LineWidth', 2)
    ylim([0 1.1]);
    ylabel("Polarity");



    dim = [0.13 0.685 .3 .3];
    time = (i_datapoint - start_point) * params.time_per_datapoint;
    str = sprintf('Time: %#.1f minutes', time/60);
    annotation('textbox', dim, 'String', str, 'FitBoxToText', 'on', 'EdgeColor', 'none');
    dim2 = [0.8 0.685 .3 .3];
    str2 = sprintf('N = %#i MTs total', mt_num(i_datapoint));
    annotation('textbox', dim2, 'String', str2, 'FitBoxToText', 'on', 'EdgeColor', 'none');

    drawnow();
    writeVideo(v, getframe(gcf));
end

close(v);
