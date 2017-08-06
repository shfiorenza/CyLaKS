clear all
n_datapoints = 10000;
length_of_microtubule = 1000;
ID = 58;

final_trajectory = zeros([n_datapoints 1]);

fileDirectory = '/home/shane/Projects/overlap_analysis/new_mt_single/%s';
fileName = 'testID.file';

data_file = fopen(sprintf(fileDirectory, fileName));
mt_array = fread(data_file, [length_of_microtubule, n_datapoints], '*int');
fclose(data_file);

for t=1:1:n_datapoints 
    for site=1:1:length_of_microtubule
       if mt_array(site, t) == ID
          final_trajectory(t, 1) = site;
          break;
       elseif site == length_of_microtubule
          final_trajectory(t, 1) = -5;
       end
    end
end

%%plot fig%%
fig1 = figure(1);
set(fig1,'Position', [50, 50, 3*480, 3*270])
plot(linspace(0, n_datapoints, n_datapoints), final_trajectory);

%%style stuff%%
%grid on
%grid minor
xlabel({'Datapoint number', sprintf('[Motor ID: %d]', ID)});
ylabel('Site on microtubule');
axis = gca;
axis.YLim = [0 length_of_microtubule];
axis.TickDir = 'out';
axis.Box = 'off';
axis.GridLineStyle = '-';
set(findall(axis, 'Type', 'Line'), 'LineWidth', 2);
