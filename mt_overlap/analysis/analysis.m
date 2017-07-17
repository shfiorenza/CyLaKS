clear all
n_timesteps = 10000;
length_of_microtubule = 1000;
n_l = 0;
n_l_temp = 0;
n_r = 0;
n_r_temp = 0;

raw_data = zeros([length_of_microtubule 2*n_timesteps]);
temp_one = zeros([length_of_microtubule 1]);
temp_two = zeros([length_of_microtubule 1]);
mt_one = zeros([length_of_microtubule 1]);
mt_two = zeros([length_of_microtubule 1]);
final_data = zeros([length_of_microtubule 2]);

data_file = fopen('/home/shane/Projects/mt_overlap/test.file');
raw_data = fread(data_file, [length_of_microtubule, 2*n_timesteps], '*int');
fclose(data_file);

raw_data((raw_data ~= 2) | (raw_data == 0)) = 0;
raw_data((raw_data == 2) | (raw_data ~= 0)) = 1;

for i=1:2:((2*n_timesteps)-1)
    temp_one(:, 1) = raw_data(:, i);
    mt_one(:, 1) = mt_one(:, 1) + double(temp_one(:, 1)./n_timesteps);
end
for i=2:2:(2*n_timesteps)
    temp_two(:, 1) = raw_data(:, i);
    mt_two(:, 1) = mt_two(:, 1) + double(temp_two(:, 1)./n_timesteps);
end

final_data(:, 1) = mt_one(:, 1);
final_data(:, 2) = mt_two(:, 1);

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
set(fig1,'Position', [50, 50, 1.5*560, 1.5*420])
plot(linspace(0, 1, length_of_microtubule), final_data);

%%style stuff%%
grid on
grid minor
xlabel('Fractional length of microtubule')
ylabel('Fraction of the time occupied')
axis = gca;
axis.YLim = [0 1];
axis.XLim = [-0.025 1.025];
axis.TickDir = 'out';
axis.Box = 'off';
axis.GridLineStyle = '-';