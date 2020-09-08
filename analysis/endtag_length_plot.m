%{
clear all
% Pseudo-constant variables
motor_speciesID = 2;
xlink_speciesID = 1;
n_datapoints = 100000;
starting_point = 50000;
active_datapoints = n_datapoints - starting_point;
fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
fileStruct = '%s_occupancy.file';
legendLabel = {'Motors', 'Crosslinkers'};

mt_lengths = [2, 4, 6, 8, 10, 14]; % in microns
xlink_concs = [0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]; % in nanomolar

n_lengths = length(mt_lengths);
n_concs = length(xlink_concs);

final_data = zeros([n_concs n_lengths]);

for i_conc=1:1:n_concs
    xlink_conc = xlink_concs(i_conc);
    for i_length=1:1:n_lengths
        n_sites = mt_lengths(i_length) * 125;
        simName = sprintf('Endtag_%#.1f_%i', xlink_conc, n_sites);
        fileName = sprintf(fileDirectory, sprintf(fileStruct, simName));
        
        motor_avg_occupancy = zeros([n_sites 1]);
        
        data_file = fopen(fileName);
        motor_raw_data = fread(data_file, [n_sites, n_datapoints], '*int');
        fclose(data_file);
        
        motor_raw_data(motor_raw_data ~= motor_speciesID) = 0;
        motor_raw_data(motor_raw_data == motor_speciesID) = 1;

        % Read in and average occupancy data over all datapoints
        for i=starting_point:1:n_datapoints
            motor_avg_occupancy(:,1) = motor_avg_occupancy(:,1) + double(motor_raw_data(:,i))./active_datapoints;
        end        
        motor_smoothed_avg = smoothdata(motor_avg_occupancy);

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
        final_data(i_conc, i_length) = endtag_length; 
    end
end
%}
fig1 = figure(1);
set(fig1,'Position', [50, 50, 2*480, 2*300])
hold on
for i=1:1:n_concs
    plot(mt_lengths, final_data(i,:),'LineWidth',2);
end

legendLabel = cellstr(num2str(xlink_concs', '%#.1f nM PRC1'));
title('Endtag length scaling for 1.5 nM kinesin-1 and PRC1 (no tethering)');
xlabel('Length of microtubule (microns)');
ylabel('Endtag length (microns)');
ylim([0 0.3]);
xlim([mt_lengths(1) mt_lengths(n_lengths)]);
legend(legendLabel, 'location', 'northwest');