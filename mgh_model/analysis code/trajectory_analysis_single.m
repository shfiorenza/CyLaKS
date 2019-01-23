clear all
n_mts = 1;
length_of_microtubule = 500;
n_datapoints = 100000;

starting_point = 00001;
n_steps = 10000;
ID = 1;

final_trajectory = zeros([n_steps 1]);

fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
fileName = 'test_motorID.file';

data_file = fopen(sprintf(fileDirectory, fileName));
mt_array = fread(data_file, [n_mts*length_of_microtubule, n_datapoints], '*int');
fclose(data_file);

% Find where motor
n_start = 1;
done = false;
for n = starting_point:1:n_datapoints
   for site = 1:1:n_mts*length_of_microtubule
       if done == false
         if mt_array(site, n) == ID
            n_start = n;
            done = true;
         end
       end
   end
end
n_end = n_start + n_steps - 1;
if(n_end > n_datapoints)
    n_end = n_datapoints;
end 

% Collect desire amount of datapoints
step_index = 1;
for n = n_start:1:n_end 
    for site = 1:1:n_mts*length_of_microtubule
            if mt_array(site, n) == ID
                % Convert raw site # to normalized index
                site_adj = mod(site, length_of_microtubule);
                final_trajectory(step_index, 1) = site_adj;
                step_index = step_index + 1; 
                break;
            % If not found, put coord that won't show up on graph
            elseif site == n_mts*length_of_microtubule
                final_trajectory(step_index, 1) = -5;
                step_index = step_index + 1;
            end
    end
end

%%plot fig%%
fig1 = figure(1);
set(fig1,'Position', [50, 50, 3*480, 3*270])
plot1 = plot(linspace(n_start, n_end, n_steps), final_trajectory, ':ob');

%%style stuff%%
plot1.MarkerSize = 2;
plot1.LineWidth = 1;
axis = gca;
axis.YLim = [0 length_of_microtubule];
axis.XLim = [n_start n_end];
axis.TickDir = 'out';
axis.Box = 'off';
axis.GridLineStyle = '-';
xlabel({'Datapoint number', sprintf('[Motor ID: %d]', ID)});
ylabel('Site on microtubule');

