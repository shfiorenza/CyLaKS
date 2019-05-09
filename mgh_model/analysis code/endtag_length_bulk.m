clear all
%n_sites = 1000;
k_on = 0.020;
% Pseudo-constant variables
motor_speciesID = 2;
xlink_speciesID = 1;
n_datapoints = 10000;
starting_point = 9000;
active_datapoints = n_datapoints - starting_point;
fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/newProc_250x_mid/%s';
fileStruct = '%s_occupancy.file';

experiment = [1.1, 1.2, 1.35, 1.55, 1.8, 2.2];

mt_lengths = [2, 4, 6, 8, 10, 14]; % in microns
runlengths = [5.0];
%runlengths = [1.2, 6.0, 12.0, 18.0, 36.0];
%motor_concs = [20,100,200
%xlink_concs = [0.5,1,2];        % in nanomolar

n_lengths = length(mt_lengths);
n_runlengths = length(runlengths);
%n_mot_concs = length(motor_concs);
%n_xl_concs = length(xlink_concs);

%final_data = zeros([n_xl_concs n_lengths]);
%final_data = zeros([n_xl_concs n_mot_concs]);
final_data = zeros([n_runlengths n_lengths]);

for i_runlength = 1:1:n_runlengths
    processivity = runlengths(i_runlength);
%for i_xl_conc=1:1:n_xl_concs
    %xlink_conc = xlink_concs(i_xl_conc);
    %for i_mot_conc=1:1:n_mot_concs   
        %motor_conc = motor_concs(i_mot_conc);
        %simName = sprintf('Endtag_%#.1fx_%i',xlink_conc, n_sites);
    for i_length=1:1:n_lengths
        n_sites = mt_lengths(i_length) * 125;  
        simName = sprintf('Endtag_%i', n_sites);
        %simName = sprintf('Endtag_%#.1f_%#.4f_%i', processivity, k_on, n_sites);
        %simName = sprintf('Endtag_%#.1f_%i', processivity, n_sites);
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
        occupancy_slope = abs(gradient(net_smoothed_avg, 0.08));
        max_slope = max(occupancy_slope);
        
        endtag_site = 0;
        past_max = false;
        for i_site=1:n_sites
            if(occupancy_slope(i_site) >= max_slope)
                past_max = true;
            end
            if(occupancy_slope(i_site) < max_slope/2 && past_max)
                endtag_site = i_site;
                break;
            end
        end
        endtag_length = endtag_site*0.008;  
        %final_data(i_xl_conc,i_length) = endtag_length;
        %final_data(i_xl_conc, i_mot_conc) = endtag_length;
        final_data(i_runlength, i_length) = endtag_length;
    end
end

fig1 = figure(1);
set(fig1,'Position', [50, 50, 2*480, 2*300])
hold on
for(i=1:1:n_runlengths)
    plot(mt_lengths, final_data(i, :),'LineWidth',2);
end
plot(mt_lengths, experiment, 'LineWidth', 2);

%legendLabel = cellstr(num2str(xlink_concs', '%#.1f nM PRC1'));
%legendLabel = cellstr(num2str(runlengths', '%#.1f um run length'));
legendLabel = {'Simulation', 'Experiment'};
%title('Endtag length scaling for 1.5 nM kinesin & 3 nM K_D for tethering');
title({'Endtag length scaling for for 5.3 micron processivity & 420 nm/s velocity'});%, ...
    %'Hydrolysis occurs 250 times slower in stalled motors'});
%title(sprintf('End-tag formation for adjusted processivity - k_{on} = %#.4f', k_on));
xlabel('Length of microtubule (microns)');
%xlabel('kinesin-1 concentration (nM)');
ylabel('Endtag length (microns)');
ylim([0 2.5]);
xlim([mt_lengths(1) mt_lengths(n_lengths)]);
legend(legendLabel, 'location', 'northwest', 'FontSize', 14);