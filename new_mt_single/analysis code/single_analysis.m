clear all
n_timesteps = 100000;
length_of_microtubule = 180;

temp_one = zeros([length_of_microtubule 1]);
mt_one = zeros([length_of_microtubule 1]);
final_data = zeros([length_of_microtubule 1]);

polarityArray = {'Plus-end on right'};

fileDirectory = '/home/shane/Projects/overlap_analysis/new_mt_overlap/%s';
fileName = 'test.file';

data_file = fopen(sprintf(fileDirectory, fileName));
raw_data = fread(data_file, [length_of_microtubule, n_timesteps], '*int');
fclose(data_file);

raw_data((raw_data ~= 2) | (raw_data == 0)) = 0;
raw_data((raw_data == 2) | (raw_data ~= 0)) = 1;

for i=1:1:n_timesteps
    temp_one(:, 1) = raw_data(:, i);
    mt_one(:, 1) = mt_one(:, 1) + double(temp_one(:, 1)./n_timesteps);
end

final_data(:, 1) = mt_one(:, 1);

%%plot fig%%
fig1 = figure(1);
set(fig1,'Position', [50, 50, 2.5*480, 2.5*300])
plot(linspace(0, 1, length_of_microtubule), final_data);

%%style stuff%%
grid on
grid minor
xlabel({'Fractional length of microtubule', sprintf('(N Sites = %d)', length_of_microtubule)});
ylabel('Fraction of the time occupied');
axis = gca;
axis.YLim = [0 1];
axis.XLim = [0 1];
axis.TickDir = 'out';
axis.Box = 'off';
axis.GridLineStyle = '-';
set(findall(axis, 'Type', 'Line'), 'LineWidth', 2);
legend(polarityArray, 'Location', 'northeastoutside');
