
clear all;
% Often-changed variables
n_sites = 1250;
simName = 'mot_15';
% Pseudo-constant variables
delta_t = 0.0005;
n_steps = 20000000;
motor_speed = 0.65;     %in um/s
smallest_time = 0.3;    %in seconds
smallest_run = motor_speed * smallest_time;  % in um
n_mts = 1;
n_datapoints = 100000;
time_per_datapoint = delta_t * n_steps / n_datapoints;
starting_point = 1;
active_datapoints = n_datapoints - starting_point;
site_size = 0.008; % in um

fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
motorFileStruct = '%s_motorID.file';
motorFileName = sprintf(fileDirectory, sprintf(motorFileStruct, simName));

motor_data_file = fopen(motorFileName);
raw_motor_data = fread(motor_data_file, [n_mts * n_sites * n_datapoints], '*int'); 
fclose(motor_data_file);
motor_data = reshape(raw_motor_data, n_sites, n_mts, n_datapoints);

starting_site = zeros([n_mts n_mts*n_sites]) - 1;
run_lengths = zeros([(n_mts * n_sites * n_datapoints) 1]);
n_runs = 1; 

for i_data = starting_point:1:n_datapoints - 1
   for i_mt = 1:1:n_mts
       motor_IDs = motor_data(:, i_mt, i_data);
       % Scan through IDs of bound motors (-1 means no motor on that site)
       for i_site = 1:1:n_sites
          motor_ID = motor_IDs(i_site);
          % If a motor is found, check if it's already been accounted for
          if motor_ID ~= -1 && motor_ID ~= 0
              % Record the motor's starting site if this is the first time
              % seeing it (-1 means it was not seen last datapoint)
              if starting_site(i_mt, motor_ID) == -1
                  starting_site(i_mt, motor_ID) = i_site;
              end         
          end
       end
       % Check one datapoint into the future to see if any motors unbound
       future_IDs = motor_data(:, i_mt, i_data + 1);
       for i_motor = 1:1:n_mts * n_sites
           motor_ID = i_motor;
           future_site = find(future_IDs == motor_ID);
           if isempty(future_site) && starting_site(i_mt, motor_ID) ~= -1
               end_site = find(motor_IDs == motor_ID);
               start_site = starting_site(i_mt, motor_ID); 
               delta = abs(end_site(1) - start_site); 
               run_length = delta * site_size;
               if run_length >= smallest_run
                   run_lengths(n_runs) = run_length;
                   n_runs = n_runs + 1;
                  
               end
               starting_site(i_mt, motor_ID) = -1;
           end  
       end
   end
end


run_lengths = run_lengths(1:n_runs - 1);

n_bins = int32(sqrt(length(run_lengths)));
hf = histfit(run_lengths, n_bins, 'exponential');

h = get(gca,'Children');
get(hf(1)); %properties of the histogram
get(hf(2)); %properties of the normal curve
%you can retrieve and plot the curve
figure;
x=get(hf(2),'XData'); 
y=get(hf(2),'YData');

max_height = max(y)

y = y / max_height;

plot(x,y);
grid on
grid minor