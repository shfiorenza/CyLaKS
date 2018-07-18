clear all

mt_lengths = [2, 4, 6, 8, 10, 14];

n_datapoints = length(mt_lengths);

final_data_00 = [n_datapoints 1];
final_data_00(1) = 0.872;
final_data_00(2) = 1.296;
final_data_00(3) = 1.24;
final_data_00(4) = 1.184;
final_data_00(5) = 1.256;
final_data_00(6) = 1.288;

final_data_02 = [n_datapoints 1];
final_data_02(1) = 0.872;
final_data_02(2) = 1.024;
final_data_02(3) = 1.368;
final_data_02(4) = 1.488;
final_data_02(5) = 1.832;
final_data_02(6) = 2.096;

final_data_04 = [n_datapoints 1];
final_data_04(1) = 0.856;
final_data_04(2) = 1.872;
final_data_04(3) = 2.664;
final_data_04(4) = 3.736;
final_data_04(5) = 4.264;
final_data_04(6) = 5.4;

plot(mt_lengths, final_data_00, 'LineWidth', 2, 'Color', 'blue');
hold on
%plot(mt_lengths, final_data_02, 'LineWidth', 2, 'Color', 'red');
%plot(mt_lengths, final_data_04, 'LineWidth', 2, 'Color', 'black');

xlabel('Length of microtubule (microns)');
ylabel('Endtag length (microns)');
ylim([0 6]);
xlim([0 15]);
legend('0.0 nM PRC1', '0.1 nM PRC1', '0.4 nM PRC1', 'location', 'northwest');

set(gcf, 'color', 'white');