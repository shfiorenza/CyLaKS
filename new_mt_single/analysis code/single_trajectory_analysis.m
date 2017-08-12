clear all
n_datapoints = 2000;
length_of_microtubule = 1000;
ID = 998;

final_trajectory = zeros([n_datapoints 1]);

fileDirectory = '/home/shane/Projects/overlap_analysis/new_mt_single/%s';
dataFileName = 'testID.file';
tubulinFileName = 'tubulin_distribution_0.file';

% Tubulin distribution along microtubule: 0 corresponds to a 
% normal tubulin site, 1 corresponds to a mutant site
tubulin_file = fopen(sprintf(fileDirectory, tubulinFileName));
tubulin_dist = fread(tubulin_file, [length_of_microtubule], '*int');
fclose(tubulin_file);

data_file = fopen(sprintf(fileDirectory, dataFileName));
mt_array = fread(data_file, [length_of_microtubule, n_datapoints], '*int');
fclose(data_file);

start = false; 
for t=1:1:n_datapoints 
    for site=1:1:length_of_microtubule
        if start == false
            if mt_array(site, t) == ID
                final_trajectory(t, 1) = site;
                n_start = t; 
                start = true;
                break;
            elseif site == length_of_microtubule
                n_datapoints = n_datapoints + 1;
            end
        elseif start == true
            if mt_array(site, t) == ID
                final_trajectory(t, 1) = site;
                break;
            elseif site == length_of_microtubule
                final_trajectory(t, 1) = -5;
            end
        end
    end
end

data_range = n_datapoints - (n_start - 1);

%%plot fig%%
fig1 = figure(1);
set(fig1,'Position', [50, 50, 3*480, 3*270])
plot1 = plot(linspace(n_start, n_datapoints, data_range), final_trajectory, ':ob')
% Plot red horizontal lines corresponding to the location of mutant sites
hold on;
for index=1:1:length_of_microtubule
    length = n_datapoints*tubulin_dist(index);
    plot([0 length], [index index], ':r', 'LineWidth', 1.0);
end 
hold off

%%style stuff%%
plot1.MarkerSize = 1.5;
plot1.LineWidth = 1;
axis = gca;
axis.YLim = [0 length_of_microtubule];
axis.XLim = [n_start n_datapoints];
axis.TickDir = 'out';
axis.Box = 'off';
axis.GridLineStyle = '-';
xlabel({'Datapoint number', sprintf('[Motor ID: %d]', ID)});
ylabel('Site on microtubule');
