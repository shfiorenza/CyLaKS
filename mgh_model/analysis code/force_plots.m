clear all;
% Often-changed variables
n_sites = 500;
simName = 'slide_scan/unteth_0.050';
% Pseudo-constant variables
n_mts = 2;
n_steps = 10000000;
n_datapoints = 10000;
delta_t = 0.000005;
starting_point = 0;
unpin_time = 10;           % time at which top MT is able to slide (in sec)
start_time = starting_point * delta_t;
end_time = n_steps * delta_t;
length = n_sites * 0.008;

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
set(fig1, 'Position', [0, 0, 1600, 800]);

sub1 = subplot(3, 1, 1);
plot(linspace(start_time, end_time, n_datapoints), xlink_data(1, :));
line([unpin_time unpin_time], ylim, 'Linestyle', '--', 'Color', 'red');
title('Xlink forces');
ylabel('Net force (pN)');
xlim([start_time end_time]);
grid on
grid minor
axis tight

subplot(3, 1, 2)
plot(linspace(start_time, end_time, n_datapoints), motor_data(1, :));
line([unpin_time unpin_time], ylim, 'Linestyle', '--', 'Color', 'red');
title('Motor (tether) forces');
ylabel('Net force (pN)');
xlim([start_time end_time]);
grid on
grid minor
axis tight

subplot(3, 1, 3)
plot(linspace(start_time, end_time, n_datapoints), total_data(1, :));
line([unpin_time unpin_time], ylim, 'Linestyle', '--', 'Color', 'red');
title('Total forces');
xlabel('Time (sec)');
ylabel('Net force (pN)');
xlim([start_time end_time]);
grid on
grid minor
axis tight

%legend({'Bottom MT', 'Top MT'});