clear all;
% Often-changed variables
n_sites = [1000,500];
%simName = '2019_11_11_slideScan/slide_scan_500_5000';
simName = 'test';
% Pseudo-constant variables
n_mts = length(n_sites);
n_steps = 6000000;
n_datapoints = 10000;
delta_t = 0.000025;
starting_point = 0001;
unpin_time = 50;           % time at which top MT is able to slide (in sec)
start_time = (starting_point/n_datapoints) * n_steps * delta_t;
end_time = n_steps * delta_t;

fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
motorFileName = '%s_motor_force.file';
xlinkFileName = '%s_xlink_force.file';
totalFileName = '%s_total_force.file';

motor_data_file = fopen(sprintf(fileDirectory, sprintf(motorFileName, simName)));
motor_data = fread(motor_data_file, [n_mts, n_datapoints], 'double');
fclose(motor_data_file);
xlink_data_file = fopen(sprintf(fileDirectory, sprintf(xlinkFileName, simName)));
xlink_data = fread(xlink_data_file, [n_mts, n_datapoints], 'double');
fclose(xlink_data_file);
total_data_file = fopen(sprintf(fileDirectory, sprintf(totalFileName, simName)));
total_data = fread(total_data_file, [n_mts, n_datapoints], 'double');
fclose(total_data_file);

fig1 = figure();
set(fig1, 'Position', [100, 0, 1200, 750]);

smoothed_xlinks = smooth(xlink_data(1,:), 100);
sub1 = subplot(3, 1, 1);
plot(linspace(start_time, end_time, n_datapoints), smoothed_xlinks, 'LineWidth', 2);
line([unpin_time unpin_time], ylim, 'Linestyle', '--', 'Color', 'red');
title('Xlink forces');
ylabel('Net force (pN)');
xlim([start_time end_time]);
grid on
grid minor
%axis tight

smoothed_motors = smooth(motor_data(1, :), 100);
subplot(3, 1, 2)
plot(linspace(start_time, end_time, n_datapoints), smoothed_motors, 'LineWidth', 2);
line([unpin_time unpin_time], ylim, 'Linestyle', '--', 'Color', 'red');
title('Motor (tether) forces');
ylabel('Net force (pN)');
xlim([start_time end_time]);
grid on
grid minor
%axis tight

smoothed_total = smooth(total_data(1,:), 100);
subplot(3, 1, 3)
plot(linspace(start_time, end_time, n_datapoints), smoothed_total, 'LineWidth', 2);
line([unpin_time unpin_time], ylim, 'Linestyle', '--', 'Color', 'red');
title('Total forces');
xlabel('Time (sec)');
ylabel('Net force (pN)');
xlim([start_time end_time]);
grid on
grid minor
%axis tight

%legend({'Bottom MT', 'Top MT'});