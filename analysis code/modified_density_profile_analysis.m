%rehauled analysis code (by Shane)
clear all;
n_timesteps = 10000;   
mt_length = 1000;
n_max_pp = 2001;      %shorthand for n_sites_max + 1; corresponds to same thing

fileDirectory = '/home/shane/Projects/overlap_analysis/%s';
fileName = 'detail2.file';

combined_mt_half_one = zeros([mt_length/2 1]);
combined_mt_half_two = zeros([mt_length/2 1]);
combined_mt_full_corrected = zeros([mt_length 1]);

mt_one_final = zeros([mt_length 1]);
mt_two_final = zeros([mt_length 1]);
mt_one_temp = zeros([mt_length 1]);
mt_two_temp = zeros([mt_length 1]);

fid1= fopen(sprintf(fileDirectory, fileName));
run_data = fread(fid1, [2*n_max_pp, n_timesteps], '*int');

for ii=1:1:n_timesteps
    % Data is written to file in terms of N_SITES_MAX, so every site beyond
    % actual length of MT is zero; need to skip them
    combined_mt_half_one(:, 1) = run_data(1:1:mt_length/2, ii);
    combined_mt_half_two(:, 1) = run_data(n_max_pp+1:1:n_max_pp + mt_length/2, ii);
    % Reverse order of first half of data, because it's written that way
    combined_mt_full_corrected(1:1:mt_length/2, 1) = combined_mt_half_one(mt_length/2:-1:1, 1);
    combined_mt_full_corrected((mt_length/2 + 1):1:mt_length, 1) = combined_mt_half_two(:, 1);
    % The data is written as one microtubule array, with the occupancy of
    % the two different MTs encoded in the integer value: 2 corresponds to
    % right-ward moving motors (on MT one), 3 corresponds to left-ward
    % moving motors (on MT two), and 4 corresponds to both (i.e. a motor is
    % bound to both microtubules at that site). This step separates the two
    mt_one_temp = combined_mt_full_corrected(:, 1);
    mt_two_temp = combined_mt_full_corrected(:, 1);
    
    mt_one_temp((mt_one_temp~=2) &(mt_one_temp~=4)) = 0;
    mt_one_temp((mt_one_temp==2) | (mt_one_temp==4)) = 1;
    mt_two_temp((mt_two_temp~=3) &(mt_two_temp~=4)) = 0;
    mt_two_temp((mt_two_temp==3) | (mt_two_temp==4)) = 1;
    % Adds this timestep to final data average
    mt_one_final(:, 1) = mt_one_final(:, 1) + double(mt_one_temp(:, 1)./n_timesteps);
    mt_two_final(:, 1) = mt_two_final(:, 1) + double(mt_two_temp(:, 1)./n_timesteps);
  
end
  
fclose(fid1);

fig1 = figure(1);
set(fig1,'Position', [50, 50, 1.5*560, 1.5*420])
plot(linspace(0, 1, mt_length), mt_one_final);
    hold on;
plot(linspace(0, 1, mt_length), mt_two_final);

ylim([0 1])
grid on
grid minor
