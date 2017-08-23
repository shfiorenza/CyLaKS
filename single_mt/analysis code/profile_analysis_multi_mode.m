clear all
n_timesteps = 60;         %in millions 
n_datapoints = 100000;	
length_of_microtubule = 180;

% Directory in which sim files are contained
fileFormat = '/home/shane/Desktop/sims/%dmil/%d/%s.file';
% Names of the different mode files
modes = {'LHLH', 'Hn', 'Hx', 'LHx', 'M', 'Ln', 'LHn', 'Lx'}; 
n_modes = numel(modes);

% Set color order so that mode colors correspond with Hui-Shun's
ColorArray = {'orange', 'cyan', 'blue', 'green', 'yellow', 'dark grey','magenta', 'red'};
RGBArray = [0.91 0.41 0.17; 0 1 1; 0 0 1; 0 0.5 0; 1 1 0; 0.3 0.3 0.3; 1 0 1; 1 0 0];

final_data = zeros([length_of_microtubule n_modes]);

% Run through all the different simulation files
for i=1:1:n_modes
    
	mt_one = zeros([length_of_microtubule 1]);
    temp_one = zeros([length_of_microtubule 1]);
    
	% Accesses simulation file that corresponds to appropriate mode
	file_name = sprintf(fileFormat, n_timesteps, length_of_microtubule, modes{i});
	data_file = fopen(file_name);
	raw_data = fread(data_file, [length_of_microtubule, 2*n_datapoints], '*int');
	fclose(data_file);
    
	% Read in data for one microtubule only
	raw_data((raw_data ~= 2) | (raw_data == 0)) = 0;
	raw_data((raw_data == 2) | (raw_data ~= 0)) = 1;
	for j=1:2:((2*n_datapoints)-1)
    	temp_one(:, 1) = raw_data(:, j);
     	mt_one(:, 1) = mt_one(:, 1) + double(temp_one(:, 1)./n_datapoints);
    end
    	
	% Transfer microtubule data to final array
	final_data(:, i) = mt_one(:, 1);
    
    clear mt_one;
    clear temp_one;
end

fig1 = figure(1);
set(0, 'DefaultAxesColorOrder', RGBArray);
plot(linspace(0, 1, length_of_microtubule), final_data);
set(fig1,'Position', [50, 50, 2.5*480, 2.5*300])

%style stuff
grid on
grid minor
xlabel({'Fractional length of microtubule', sprintf('(N Sites = %d)', length_of_microtubule)});
ylabel('Fraction of the time occupied');
axis = gca;
axis.YLim = [0 1];
axis.XLim = [0 1];
axis.TickDir = 'out';
axis.Box = 'off';
axis.GridLineStyle = '-';
set(findall(axis, 'Type', 'Line'), 'LineWidth', 2);
legend(modes, 'Location', 'northeastoutside');



