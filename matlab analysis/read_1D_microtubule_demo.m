clear all;
n_timesteps = 1000;   
length_of_mt = 1000;

fileDirectory = '/home/shane/Projects/overlap_analysis/%s';
fileName = 'testdetail.file';

raw_mt_data = zeros([length_of_mt 1]);
corrected_mt_data = zeros([length_of_mt 1]);

mt_one = zeros([length_of_mt 1]);
mt_two = zeros([length_of_mt 1]);
mt_one_temp = zeros([length_of_mt 1]);
mt_two_temp = zeros([length_of_mt 1]);

fid1= fopen(sprintf(fileDirectory, fileName));
run_data = fread(fid1, [length_of_mt, n_timesteps], '*int');

for ii=1:1:n_timesteps
    raw_mt_data(:, 1) = run_data(:, ii);
        
    % Reverses order of first half of data, because apparently that's
    % how it's written to file
    corrected_mt_data(1:1:length_of_mt/2, 1) = raw_mt_data(length_of_mt/2:-1:1, 1);
    corrected_mt_data(length_of_mt/2+1:1:length_of_mt, 1) = raw_mt_data(length_of_mt/2+1:1:length_of_mt, 1);
    
    mt_one_temp = corrected_mt_data(:, 1);
    mt_two_temp = corrected_mt_data(:, 1);
    
    % The data is written as one microtubule array, with the occupancy of
    % the two different MTs encoded in the integer value: 2 corresponds to
    % right-ward moving motors (on MT one), 3 corresponds to left-ward
    % moving motors (on MT two), and 4 corresponds to both (i.e. a motor is
    % bound to both microtubules at that site). This step separates them
    mt_one_temp((mt_one_temp~=2) &(mt_one_temp~=4)) = 0;
    mt_one_temp((mt_one_temp==2) | (mt_one_temp==4)) = 1;
    mt_two_temp((mt_two_temp~=3) &(mt_two_temp~=4)) = 0;
    mt_two_temp((mt_two_temp==3) | (mt_two_temp==4)) = 1;
    
    % Adds this timestep to final data average
    mt_one(:, 1) = mt_one(:, 1) + double(mt_one_temp(:, 1)./n_timesteps);
    mt_two(:, 1) = mt_two(:, 1) + double(mt_two_temp(:, 1)./n_timesteps);
  
end
  
fclose(fid1);

fig1 = figure(1);
set(fig1,'Position', [50, 50, 1.5*560, 1.5*420])

plot(linspace(0, 1, length_of_mt), mt_one); %right
    hold on;
plot(linspace(0, 1, length_of_mt), mt_two); %left
ylim([0 1])