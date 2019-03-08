clear all
% Pseudo-constant variables
motor_speciesID = 2;
xlink_speciesID = 1;
n_datapoints = 10000;
starting_point = 5000;
active_datapoints = n_datapoints - starting_point;
fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/scan_endtag_3nM/%s';
fileStruct = '%s_occupancy.file';

mt_lengths = [2, 4, 6, 8, 10, 14]; % in microns
xlink_concs = [0.0, 0.1, 0.4]; % in nanomolar

n_lengths = length(mt_lengths);
n_concs = length(xlink_concs);

final_data = zeros([n_concs n_lengths]);

%sim_names = ["Endtag_00.0x_%i", "Endtag_0.1.0x_%i", "Endtag_0.4.0x_%i"];

for i_conc=1:1:n_concs
    xlink_conc = xlink_concs(i_conc);
    for i_length=1:1:n_lengths
        n_sites = mt_lengths(i_length) * 125;
        %simName = sprintf(sim_names(i_conc), n_sites);
        simName = sprintf('Endtag_%#.1fx_%i',xlink_conc, n_sites);
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
        for i=1:n_sites
            site_occupancy = net_smoothed_avg(i);
            if(site_occupancy >= 0.5)
                endtag_site = i;
                break;
            end
        end
        endtag_length = endtag_site*0.008;
        
        
        final_data(i_conc,i_length) = endtag_length;
    end
end

fig1 = figure(1);
set(fig1,'Position', [50, 50, 2*480, 2*300])
hold on
for i=1:1:n_concs
    plot(mt_lengths, final_data(i,:),'LineWidth',2);
end

legendLabel = cellstr(num2str(xlink_concs', '%#.1f nM PRC1'));
title('Endtag length scaling for 1.5 nM kinesin & 3 nM K_D for tethering');
xlabel('Length of microtubule (microns)');
ylabel('Endtag length (microns)');
ylim([0 0.3]);
xlim([mt_lengths(1) mt_lengths(n_lengths)]);
legend(legendLabel, 'location', 'best');