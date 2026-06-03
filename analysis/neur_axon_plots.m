sim_name = 'test_long';

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


mt_num = filament_num_raw; 
mt_len = NaN(params.n_datapoints, max(mt_num));
mt_pos = NaN(params.n_datapoints, max(mt_num), 2, params.n_dims);

i_entry = 1;
j_entry = 1;
for i_datapoint = 1 : 1 : params.n_datapoints
    for i_mt = 1 : mt_num(i_datapoint)
        mt_len(i_datapoint, i_mt) = filament_lengths_raw(i_entry);
        i_entry = i_entry + 1;
        for i_end = 1 : 1 : 2
            for i_dim = 1 : 1 : params.n_dims
                mt_pos(i_datapoint, i_mt, i_end, i_dim) = filament_pos_raw(j_entry);
                j_entry = j_entry + 1;
            end
        end
    end
end

%pos = zeros(params.n_datapoints, 10);
minus_pos = NaN(params.n_datapoints, max(mt_num));
for i_datapoint = 1 : 1 : params.n_datapoints
   for i_mt = 1 : mt_num(i_datapoint)
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


fig = figure()

plot(linspace(0, params.t_run / 60, params.n_datapoints), mt_num, 'LineWidth', 3)
fontname("arial");
fontsize(14, "points");
xlabel("Time (min)", "FontSize", 18);
ylabel("Microtubule count", "FontSize", 18);
ylim([0 700]);

return
plot(pos(:, 1));
return
hold on
for i_mt = 200 : 205
    plot(pos(:, i_mt));
end

%min_x = min(mt_pos(:, :, :, 1), [], "all");
%max_x = max(mt_pos(:, :, :, 1), [], "all");
%span = max_x - min_x;
%n_sites = span / 8.2;