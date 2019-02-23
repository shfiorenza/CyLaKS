clear all
% Often-changed variables
n_sites = 1750;
xlink_conc = 0.6;
simName = 'scan_endtag_300nM/Endtag_0.4.0x_1750';
%simName = sprintf('Endtag_%#.1f_%i', xlink_conc, n_sites);
% Pseudo-constant variables
motor_speciesID = 2;
xlink_speciesID = 1;
n_datapoints = 10000;
starting_point = 50000;
active_datapoints = n_datapoints - starting_point;

%fileDirectory = '/media/shane/Shane''s External HDD (1 TB)/Parameter Scan 1/%s';
%fileDirectory = '/home/shane/Desktop/pseudo_crackpot/%s';
fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
fileStruct = '%s_occupancy.file';
legendLabel = {'Motors', 'Crosslinkers'};

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

% Scan through smoothed occupancy data to find approx end-tag position
background = motor_smoothed_avg(n_sites);
max_occu = background;
for i=1:1:n_sites
    site_occupancy = motor_smoothed_avg(i);
    if(site_occupancy > max_occu)
        max_occu = site_occupancy;
    end
end
endtag_site = n_sites;
for i=n_sites:-1:1
    site_occupancy = motor_smoothed_avg(i);
    if(site_occupancy > max_occu / 2)
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
plot(linspace(0, n_sites*0.008, n_sites), xlink_smoothed_avg)
% Put vertical red line where endtag starting position is
plot([endtag_length endtag_length], [0 1], ':r', 'LineWidth', 0.1);

%%style stuff%%

title({sprintf('Microtubule of length %g microns', n_sites*0.008), ...
    sprintf('Endtag length: %g microns', endtag_length), ...
    sprintf('Crosslinker concentration: %g nM', xlink_conc)});
xlabel({'Distance along microtubule relative to plus-end (microns)'});
ylabel('Fraction of the time occupied');
xlim([0 100*0.008]);
ylim([0 1]);
grid on
grid minor
axis = gca;
axis.TickDir = 'out';
axis.Box = 'off';
axis.GridLineStyle = '-';
set(findall(axis, 'Type', 'Line'), 'LineWidth', 2);
legend(legendLabel, 'Location', 'northeast');