clear all
n_datapoints = 2000;
n_mts = 2;
length_of_microtubule = 1000;
ID = 358;

final_trajectory = zeros([n_datapoints 1]);

fileDirectory = '/home/shane/Projects/overlap_analysis/new_mt_overlap/%s';
fileName = 'testID.file';

data_file = fopen(sprintf(fileDirectory, fileName));
mt_array = fread(data_file, [n_mts*length_of_microtubule, n_datapoints], '*int');
fclose(data_file);

start = false; 
for t=1:1:n_datapoints 
    for site=1:1:n_mts*length_of_microtubule
        % Only show data AFTER the motor's first bind event
        if start == false
            if mt_array(site, t) == ID
                % Convert raw site % to normalized index
                site_adj = mod(site, length_of_microtubule);
                final_trajectory(t, 1) = site_adj;
                n_start = t; 
                start = true;
                break;
            elseif site == n_mts*length_of_microtubule
                n_datapoints = n_datapoints + 1;
            end
        elseif start == true
            if mt_array(site, t) == ID
                % Convert raw site # to normalized index
                site_adj = mod(site, length_of_microtubule);
                final_trajectory(t, 1) = site_adj;
                break;
            % If not found, put coord that won't show up on graph
            elseif site == n_mts*length_of_microtubule
                final_trajectory(t, 1) = -5;
            end
        end
    end
end

data_range = n_datapoints - (n_start - 1);
%%plot fig%%
fig1 = figure(1);
set(fig1,'Position', [50, 50, 3*480, 3*270])
plot1 = plot(linspace(n_start, n_datapoints, data_range), final_trajectory, ':ob');

%%style stuff%%
plot1.MarkerSize = 2;
plot1.LineWidth = 1;
axis = gca;
axis.YLim = [0 length_of_microtubule];
axis.XLim = [n_start n_datapoints];
axis.TickDir = 'out';
axis.Box = 'off';
axis.GridLineStyle = '-';
xlabel({'Datapoint number', sprintf('[Motor ID: %d]', ID)});
ylabel('Site on microtubule');

