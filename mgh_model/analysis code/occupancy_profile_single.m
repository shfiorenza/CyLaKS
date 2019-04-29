clear all
% Often-changed variables
n_sites = 1750;
processivity = 6.0;
xlink_conc = 2.0;
motor_conc = 20;
%simName = sprintf('endtag_scan/Endtag_%#.1fx_%#.1fm_%i', xlink_conc, motor_conc, n_sites);
simName = sprintf('output_end/Endtag_%#.1f_%i', processivity, n_sites);
%simName = 'Endtag_1.2_1000b';
% Pseudo-constant variables
motor_speciesID = 2;
xlink_speciesID = 1;
n_datapoints = 10000;
starting_point = 5000;
active_datapoints = n_datapoints - starting_point;

%fileDirectory = '/media/shane/Shane''s External HDD (1 TB)/Parameter Scan 1/%s';
%fileDirectory = '/home/shane/Desktop/pseudo_crackpot/%s';
fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
fileStruct = '%s_occupancy.file';
legendLabel = {'Motors', 'Crosslinkers', 'Combined'};

fileName = sprintf(fileDirectory, sprintf(fileStruct, simName));
data_file = fopen(fileName);
motor_raw_data = fread(data_file, [n_sites, n_datapoints], '*int');
xlink_raw_data = motor_raw_data;
fclose(data_file);

motor_avg_occupancy = zeros([n_sites 1]);
xlink_avg_occupancy = zeros([n_sites 1]);

motor_raw_data(motor_raw_data ~= motor_speciesID) = 0;
motor_raw_data(motor_raw_data == motor_speciesID) = 1;

xlink_raw_data(xlink_raw_data ~= xlink_speciesID) = 0;
xlink_raw_data(xlink_raw_data == xlink_speciesID) = 1;

% Read in and average occupancy data over all datapoints
for i=starting_point:1:n_datapoints
    motor_avg_occupancy(:,1) = motor_avg_occupancy(:,1) + double(motor_raw_data(:,i))./active_datapoints;
    xlink_avg_occupancy(:,1) = xlink_avg_occupancy(:,1) + double(xlink_raw_data(:,i))./active_datapoints;
end

motor_smoothed_avg = smoothdata(motor_avg_occupancy);
xlink_smoothed_avg = smoothdata(xlink_avg_occupancy);
net_smoothed_avg = motor_smoothed_avg + xlink_smoothed_avg;

endtag_site = 0;
curve = gradient(abs(gradient(motor_smoothed_avg,0.008)),0.5);
endtag_curve = max(curve);

for i=1:n_sites
    %{
    if(curve(i) == endtag_curve)
        endtag_site = i;
        break;
    end
    %}
    site_occupancy = net_smoothed_avg(i);
    if(site_occupancy < 0.3)
        endtag_site = i;
        break;
    end
end
endtag_length = endtag_site*0.008;

%%plot fig%%

fig1 = figure(1);
set(fig1,'Position', [50, 50, 2*480, 2*300])
plot(linspace(0, n_sites*0.008, n_sites), motor_smoothed_avg);

hold on
%plot(linspace(0, n_sites*0.008, n_sites), curve);
%plot(linspace(0, n_sites*0.008, n_sites), xlink_smoothed_avg);
%plot(linspace(0, n_sites*0.008, n_sites), net_smoothed_avg);

% Put vertical red line where endtag starting position is
plot([endtag_length endtag_length], [0 1], ':r', 'LineWidth', 0.1);

%%style stuff%%


%title({sprintf('%g micron-long microtubule', n_sites*0.008), ...
    %sprintf('%#.1f nM PRC1 and %#.1f nM kinesin-1', xlink_conc, motor_conc), ...
title({sprintf('Processivity: %#.1f microns', processivity), ...
    sprintf('Endtag length: %g microns', endtag_length)});

xlabel({'Distance along microtubule relative to plus-end (microns)'});
ylabel('Fraction of the time occupied');
%xlim([0 0.5]);
ylim([0 1]);
grid on
grid minor
axis = gca;
axis.TickDir = 'out';
axis.Box = 'off';
axis.GridLineStyle = '-';
set(findall(axis, 'Type', 'Line'), 'LineWidth', 2);
%legend(legendLabel, 'Location', 'northeast');