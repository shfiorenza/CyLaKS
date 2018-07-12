clear all

mt_lengths = [2, 4, 6, 8, 10, 14];

n_datapoints = length(mt_lengths);

final_data_00 = [n_datapoints 1];
final_data_00(1) = 1.312;
final_data_00(2) = 1.968;
final_data_00(3) = 0.008;
final_data_00(4) = 0.016;
final_data_00(5) = 0.008;
final_data_00(6) = 0.008;

final_data_02 = [n_datapoints 1];
final_data_02(1) = 0.072;
final_data_02(2) = 0.12;
final_data_02(3) = 0.144;
final_data_02(4) = 0.152;
final_data_02(5) = 0.144;
final_data_02(6) = 0.16;

final_data_04 = [n_datapoints 1];
final_data_04(1) = 0.136;
final_data_04(2) = 0.192;
final_data_04(3) = 0.244;
final_data_04(4) = 0.232;
final_data_04(5) = 0.24;
final_data_04(6) = 0.216;

plot(mt_lengths, final_data_00, 'LineWidth', 2, 'Color', 'blue');
hold on
%plot(mt_lengths, final_data_02, 'LineWidth', 2, 'Color', 'red');
%plot(mt_lengths, final_data_04, 'LineWidth', 2, 'Color', 'black');

xlabel('Length of microtubule (microns)');
ylabel('Endtag length (microns)');
ylim([0 0.5]);
xlim([0 15]);
legend('0.0 nM PRC1', '0.2 nM PRC1', '0.4 nM PRC1');

set(gcf, 'color', 'white');