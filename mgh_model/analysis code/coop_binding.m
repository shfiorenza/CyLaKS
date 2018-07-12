clear all;

range_one = [-9, -8, -7, -6, -5];
affinity_one = [1.5, 1.5, 1.5, 1.5, 1.5];

range_two = [-5, -4, -3];
affinity_two = [2, 2, 2];

range_three = [-3, -2, -1];
affinity_three = [2.5, 2.5, 2.5];

range_four = [-1, 0];
affinity_four = [3, 3];

range_five = [0, 1];
affinity_five = [3.5, 3.5];

range_six = [1, 2, 3];
affinity_six = [3, 3, 3,];

range_seven = [3, 4, 5];
affinity_seven = [2.5, 2.5, 2.5];

range_eight = [5, 6, 7];
affinity_eight = [2, 2, 2];

range_nine = [7, 8, 9];
affinity_nine = [1.5, 1.5, 1.5];

range = [-9, 9];
baseline = [1, 1];

center = [0 0];
line = [0 6];

hold on
plot(range, baseline, '--black', 'LineWidth', 1);
plot(center, line, '--red', 'LineWidth', 1)
plot(range_one, affinity_one, 'Color', 'black', 'LineWidth', 2);
plot(range_two, affinity_two, 'Color', 'black', 'LineWidth', 2);
plot(range_three, affinity_three, 'Color', 'black', 'LineWidth', 2);
plot(range_four, affinity_four, 'Color', 'black', 'LineWidth', 2);
plot(range_five, affinity_five, 'Color', 'black', 'LineWidth', 2);
plot(range_six, affinity_six, 'Color', 'black', 'LineWidth', 2);
plot(range_seven, affinity_seven, 'Color', 'black', 'LineWidth', 2);
plot(range_eight, affinity_eight, 'Color', 'black', 'LineWidth', 2);
plot(range_nine, affinity_nine, 'Color', 'black', 'LineWidth', 2);

ylim([0 4]);
xlim([-9 9]);
ylabel('Binding rate compared to baseline');
xlabel('Distance from bound kinesin (microns)');