clear all;
n_datapoints = 100000;
length_of_microtubule = 1000;

polarityArray = {'Plus-end on right'};

fileDirectory = '/home/shane/Projects/overlap_analysis/new_mt_single/%s';
dataFileName = 'test.file';
tubulinFileName = 'tubulin_distribution_0.file';

% Tubulin distribution along microtubule: 0 corresponds to a 
% normal tubulin site, 1 corresponds to a mutant site
tubulin_file = fopen(sprintf(fileDirectory, tubulinFileName));
tubulin_dist = fread(tubulin_file, [length_of_microtubule, 1], '*int');
fclose(tubulin_file);

n_normal = 0;
n_mutant = 0;
% Iterate through tubulin array and count up the different species
for i=1:1:length_of_microtubule
   tubulin_type = tubulin_dist(i);
   if i ~= 1 && i ~= length_of_microtubule
       if tubulin_type == 0
           n_normal = n_normal + 1;
       elseif tubulin_type == 1
           n_mutant = n_mutant + 1;
       else
           disp('what')
       end
   end
end
fprintf('%i normal, %i mutant', n_normal, n_mutant)

% Motor occupancy data: 2 corresponds to a single-headed motor
% on the site; 0 corresponds to an empty site
data_file = fopen(sprintf(fileDirectory, dataFileName));
raw_data = fread(data_file, [length_of_microtubule, n_datapoints], '*int');
fclose(data_file);

% Convert kinesin ID to binary occupancy value
raw_data((raw_data ~= 2) | (raw_data == 0)) = 0;
raw_data((raw_data == 2) | (raw_data ~= 0)) = 1;

temp_one = zeros([length_of_microtubule 1]);
final_mt = zeros([length_of_microtubule 1]);

% Iterate through all data points
for i=1:1:n_datapoints
    % Copy microtubule profile from ith data point
    temp_one(:, 1) = raw_data(:, i);
    % Average the occupancy of entire MT lattice over all data points
    final_mt(:, 1) = final_mt(:, 1) + double(temp_one(:, 1)./n_datapoints);
end

% Set figure position and size
fig1 = figure(1);
set(fig1,'Position', [50, 50, 2.5*480, 2.5*300])
% Plot final microtubule with averaged-out occupancy 
plot1 = plot(linspace(0, 1, length_of_microtubule), final_mt,'-b', 'LineWidth', 1.5);
% Plot red vertical lines corresponding to the location of mutant sites
hold on;
for index=1:1:length_of_microtubule
    
    site = index/length_of_microtubule;
    
    plot([site site], [0 tubulin_dist(index)], ':r', 'LineWidth', 1.0);
end 
hold off

%%style stuff%%
grid on
axis = gca;
axis.YLim = [0 1];
axis.XLim = [0 1];
axis.TickDir = 'out';
axis.Box = 'off';
axis.GridLineStyle = '-';
xlabel({'Fractional length of microtubule', sprintf('(%i sites: %i normal, %i mutant, 2 boundary)', length_of_microtubule, n_normal, n_mutant)});
ylabel('Fraction of the time occupied');
legend(polarityArray, 'Location', 'northeastoutside');

