clear all
n_datapoints = 6;

final_data_short = [n_datapoints 1];
final_data_short(1, :) = 0.948;
final_data_short(2, :) = 1.296;
final_data_short(3, :) = 1.368;
final_data_short(4, :) = 1.384;
final_data_short(5, :) = 1.384;
final_data_short(6, :) = 1.358;

final_data_long = [n_datapoints 1];

final_data_long(1, :) = 1.4;
final_data_long(2, :) = 1.824;
final_data_long(3, :) = 1.944;
final_data_long(4, :) = 1.992;
final_data_long(5, :) = 1.976;
final_data_long(6, :) = 1.928;

plot(linspace(2, 14, 6), final_data_short);
hold on
%plot(linspace(2, 14, 6), final_data_long);

xlabel('Length of microtubule (microns)');
ylabel('Length of endtag (microns)');