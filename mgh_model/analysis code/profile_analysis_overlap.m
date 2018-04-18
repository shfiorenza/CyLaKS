clear all
n_datapoints = 100000;
length_of_microtubule = 250;
raw_overlap_length = 125;
overlap_length = abs(raw_overlap_length);
delta = length_of_microtubule - overlap_length;
if(overlap_length == 0)
    delta = 0;
end
speciesID = 1;

n_l = 0;
n_l_temp = 0;
n_r = 0;
n_r_temp = 0;
temp_one = zeros([length_of_microtubule 1]);
temp_two = zeros([length_of_microtubule 1]);
mt_one = zeros([length_of_microtubule 1]);
mt_two = zeros([length_of_microtubule 1]);
final_data = zeros([(length_of_microtubule + delta) 2]);

polarityArray = {'Plus-end on right', 'Plus end on left'};

fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
fileName = 'presXL5_occupancy.file';

data_file = fopen(sprintf(fileDirectory, fileName));
raw_data = fread(data_file, [length_of_microtubule, 2*n_datapoints], '*int');
fclose(data_file);

raw_data((raw_data ~= speciesID) | (raw_data == 0)) = 0;
raw_data((raw_data == speciesID) | (raw_data ~= 0)) = 1;

% Copy MT data from file
for i=1:2:((2*n_datapoints)-1)
    temp_one(:, 1) = raw_data(:, i);
    mt_one(:, 1) = mt_one(:, 1) + double(temp_one(:, 1)./n_datapoints);
end
for i=2:2:(2*n_datapoints)
    temp_two(:, 1) = raw_data(:, i);
    mt_two(:, 1) = mt_two(:, 1) + double(temp_two(:, 1)./n_datapoints);
end

if(raw_overlap_length > 0)
    final_data(1:1:length_of_microtubule, 1) = mt_one(:, 1);
    final_data(delta + 1:1:(length_of_microtubule + delta), 2) = mt_two(:, 1);
    else if(raw_overlap_length < 0)
        final_data(delta + 1:1:(length_of_microtubule + delta), 1) = mt_one(:, 1);
        final_data(1:1:length_of_microtubule, 2) = mt_two(:, 1);
        else
        final_data(:, 1) = mt_one(:, 1);
        final_data(:, 2) = mt_two(:, 1);
        end;
    end
  
for i = 1:1:length_of_microtubule
    n_r_temp = final_data(i, 1);
    n_l_temp = final_data(i, 2);
    n_r = n_r + n_r_temp;
    n_l = n_l + n_l_temp;
end

%%calculate and display avg occupancy%%
n_avg = n_r + n_l;
disp(n_avg);

%%plot fig%%
fig1 = figure(1);
set(fig1,'Position', [50, 50, 2.5*480, 2.5*300])
plot(linspace(0, length_of_microtubule + delta, length_of_microtubule + delta), ...
    final_data, 'LineWidth', 1.5);

%%style stuff%%
grid on
grid minor
axis = gca;
axis.YLim = [0 1];
axis.XLim = [0 length_of_microtubule + delta];
axis.TickDir = 'out';
axis.Box = 'off';
axis.GridLineStyle = '-';
xlabel({'Coordinate', sprintf('(%d sites on each MT)', length_of_microtubule)});
ylabel('Fraction of the time occupied');
legend(polarityArray, 'Location', 'northeastoutside');
