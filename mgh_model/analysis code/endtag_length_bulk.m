clear variables;
off_ratio = 10;
k_hydrolyze = 90;
jam_ratio = 900;
% Pseudo-constant variables
mt_lengths = [2, 4, 6, 8, 10, 14]; % in microns
experiment_endtags = [1.1, 1.2, 1.35, 2.55, 1.8, 2.2];
exp_endtags = [1.15,1.25,1.45,1.6,1.8,2.15]; 
exp_errors = [0.1,0.1,0.2,0.3,0.25,0.4];
motor_speciesID = 2;
xlink_speciesID = 1;
n_datapoints = 10000;
starting_point = 5000;
active_datapoints = n_datapoints - starting_point;
fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
fileStruct = '%s_occupancy.file';


runlengths = [5.0];
%runlengths = [1.2, 6.0, 12.0, 18.0, 36.0];
n_lengths = length(mt_lengths);
n_runlengths = length(runlengths);
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
        %simName = sprintf('Endtag_%i', n_sites);     
        if(off_ratio == 1)
            simName = sprintf('scan_output/k_hydrolyze_%i/jam_ratio_%i/endtag_%i_%i_%i', ...
                k_hydrolyze, jam_ratio, k_hydrolyze, jam_ratio, n_sites);
         %  simName = sprintf('endtag_%i_%i_%i_%i_low_kOn', off_ratio, k_hydrolyze, jam_ratio, n_sites); 
        else
            simName = sprintf('scan_output_b/off_ratio_%i/k_hydrolyze_%i/jam_ratio_%i/endtag_%i_%i_%i_%i', ...
               off_ratio, k_hydrolyze, jam_ratio, off_ratio, k_hydrolyze, jam_ratio, n_sites);
            %simName = sprintf('endtag_%i_%i_%i_%i', off_ratio, k_hydrolyze, jam_ratio, n_sites);
        end
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
%plot(mt_lengths, experiment_endtags, 'LineWidth', 2);
errorbar(mt_lengths, exp_endtags, exp_errors, 'LineWidth', 2);

%legendLabel = cellstr(num2str(xlink_concs', '%#.1f nM PRC1'));
%legendLabel = cellstr(num2str(runlengths', '%#.1f um run length'));
legendLabel = {'Simulation', 'Experiment'};
%title('Endtag length scaling for 1.5 nM kinesin & 3 nM K_D for tethering');
title(sprintf('off ratio = %g, k hydrolyze = %g, jam ratio = %g', off_ratio, k_hydrolyze, jam_ratio));
%title(sprintf('End-tag formation for adjusted processivity - k_{on} = %#.4f', k_on));
xlabel('Length of microtubule (microns)');
%xlabel('kinesin-1 concentration (nM)');
ylabel('Endtag length (microns)');
%ylim([0 2.5]);
%xlim([mt_lengths(1) mt_lengths(n_lengths)]);
legend(legendLabel, 'location', 'northwest', 'FontSize', 14);