clear all
n_sites = 1000;
% Pseudo-constant variables
motor_speciesID = 2;
xlink_speciesID = 1;
n_datapoints = 10000;
starting_point = 5000;
active_datapoints = n_datapoints - starting_point;
fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
fileStruct = '%s_occupancy.file';

%mt_lengths = [2, 4, 6, 8, 10, 14]; % in microns
motor_concs = [20,100,200];
xlink_concs = [0.5,1,2];        % in nanomolar

%n_lengths = length(mt_lengths);
n_mot_concs = length(motor_concs);
n_xl_concs = length(xlink_concs);

%final_data = zeros([n_xl_concs n_lengths]);
final_data = zeros([n_xl_concs n_mot_concs]);

%sim_names = ["Endtag_00.0x_%i", "Endtag_0.1.0x_%i", "Endtag_0.4.0x_%i"];

for i_xl_conc=1:1:n_xl_concs
    xlink_conc = xlink_concs(i_xl_conc);
    %for i_length=1:1:n_lengths
    for i_mot_conc=1:1:n_mot_concs   
        %n_sites = mt_lengths(i_length) * 125;
        motor_conc = motor_concs(i_mot_conc);
        %simName = sprintf('Endtag_%#.1fx_%i',xlink_conc, n_sites);
        simName = sprintf('endtag_scan/Endtag_%#.1fx_%#.1fm_1000', xlink_conc, motor_conc);
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
            site_occupancy = xlink_smoothed_avg(i);
            if(site_occupancy < 0.05)
                endtag_site = i;
                break;
            end
        end
        endtag_length = endtag_site*0.008;
        
        %final_data(i_xl_conc,i_length) = endtag_length;
        final_data(i_xl_conc, i_mot_conc) = endtag_length;
    end
end

fig1 = figure(1);
set(fig1,'Position', [50, 50, 2*480, 2*300])
hold on
for(i=1:1:n_xl_concs)
    plot(motor_concs, final_data(i, :),'LineWidth',2);
end

legendLabel = cellstr(num2str(xlink_concs', '%#.1f nM PRC1'));
%title('Endtag length scaling for 1.5 nM kinesin & 3 nM K_D for tethering');
title('Endtag length scaling for 8-micron MT');
%xlabel('Length of microtubule (microns)');
xlabel('kinesin-1 concentration (nM)');
ylabel('Endtag length (microns)');
ylim([0 0.3]);
%xlim([mt_lengths(1) mt_lengths(n_lengths)]);
legend(legendLabel, 'location', 'best');