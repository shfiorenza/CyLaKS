
clear variables;
xlink_concs = [0, 2, 4]; 
name_structs = ["Endtag_0_%i", "Endtag_%i_02", "Endtag_%i"];
mt_lengths = [2, 4, 6, 8, 10, 14];     % in microns
% experimental parameters
exp_mt_lengths = [2.4, 4.0, 5.6, 7.2, 8.8, 10.3, 13.5];
exp_endtags_0 = [1.2, 1.3, 1.4, 1.6, 1.65, 1.75, 2.25]; 
exp_errs_0_y = [0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.2];
exp_endtags_2 = [1.4, 1.6, 1.7, 2.0, 2.1, 2.4, -1];
exp_errs_2_y = [0.1, 0.6, 0.2, 0.3, 0.2, 0.3, 0];
exp_endtags_4 = [1.6, 2.1, 2.4, 3.1, 3.6, 4.0, 5.3];
exp_errs_4_y = [0.2, 0.6, 0.25, 0.4, 0.25, 1.0, 0.25 ];
exp_endtags = [exp_endtags_0.', exp_endtags_2.', exp_endtags_4.'].';
exp_errs_x = [0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8];
exp_errs_y = [exp_errs_0_y.', exp_errs_2_y.', exp_errs_4_y.'].';
exp_line_x = [2 14.5];
exp_line_y = [[1.15 2.23]; [1.37 2.88]; [1.67 5.05]];
motor_speciesID = 2;
xlink_speciesID = 1;
n_datapoints = 10000;
starting_point = 5000;
active_datapoints = n_datapoints - starting_point;
fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/old_output/%s';
fileStruct = '%s_occupancy.file';

n_concs = length(xlink_concs);
n_lengths = length(mt_lengths);
final_data = zeros([n_concs n_lengths]);

for i_conc=1:1:n_concs
%    xlink_conc = xlink_concs(i_conc);
    for i_length=1:1:n_lengths
        n_sites = mt_lengths(i_length) * 125;  
        simName = sprintf(name_structs(i_conc), n_sites);       
        fileName = sprintf(fileDirectory, sprintf(fileStruct, simName));
        data_file = fopen(fileName);
        motor_raw_data = fread(data_file, [n_sites, n_datapoints], '*int');
        xlink_raw_data = motor_raw_data;
        fclose(data_file);
        
        motor_raw_data(motor_raw_data ~= motor_speciesID) = 0;
        motor_raw_data(motor_raw_data == motor_speciesID) = 1;  
        xlink_raw_data(xlink_raw_data ~= xlink_speciesID) = 0;
        xlink_raw_data(xlink_raw_data == xlink_speciesID) = 1;
        
        motor_avg_occupancy = zeros([n_sites 1]);
        xlink_avg_occupancy = zeros([n_sites 1]); 
        
        % Read in and average occupancy data over all datapoints
        for i=starting_point:1:n_datapoints
            motor_avg_occupancy(:,1) = motor_avg_occupancy(:,1) + double(motor_raw_data(:,i))./active_datapoints;
            xlink_avg_occupancy(:,1) = xlink_avg_occupancy(:,1) + double(xlink_raw_data(:,i))./active_datapoints;
        end

        smooth_window = n_sites / 10;   
        motor_occupancy = smoothdata(motor_avg_occupancy, 'movmean', smooth_window);
        xlink_occupancy = smoothdata(xlink_avg_occupancy, 'movmean', smooth_window);
        net_occupancy = motor_occupancy + xlink_occupancy;
        occupancy_slope = smoothdata(gradient(net_occupancy, 0.008),'movmean', smooth_window);
        occupancy_accel = smoothdata(gradient(occupancy_slope, 0.008), 'movmean', smooth_window);  
        max_occupancy = max(net_occupancy);
        min_slope = min(occupancy_slope);
  
        past_threshold = false;
        i_threshold = 0;   
        endtag_site = 0;
        for i_site=1:n_sites
            if(~past_threshold && net_occupancy(i_site) < 0.5*max_occupancy)
                past_threshold = true;
                i_threshold = i_site; 
            end
            if(past_threshold && occupancy_slope(i_site) > (min_slope / 2))
                endtag_site = i_site;
                break;
            end
        end
        if endtag_site == 0
            [max_accel, i_peak] = max(occupancy_accel(i_threshold:n_sites));
            endtag_site = i_threshold + i_peak;
        end
        endtag_length = endtag_site*0.008;  
        
        final_data(i_conc, i_length) = endtag_length;
    end
end
%}

set(0, 'DefaultLineLineWidth', 1.25);
set(0,'DefaultLegendAutoUpdate','off');
fig1 = figure(1);
set(fig1,'Position', [50, 50, 2*480, 2*300])
hold on

colors = ['b'; 'r'; 'k'];
err_config = [{'^b', 'b', 'b'}; {'^r', 'r', 'r'}; {'^k', 'k', 'k'}];
for(i=1:1:n_concs)
    plot(mt_lengths, final_data(i, :),'*','LineWidth',2, 'MarkerSize', 15, 'Color', colors(i,:));
end

%legendLabel = cellstr(num2str(xlink_concs' / 10, '%#.1f nM PRC1'));
legendLabel = {'0.0 nM PRC1', '0.1 nM PRC1', '0.4 nM PRC1'};
legend(legendLabel, 'location', 'northwest', 'FontSize', 14);

for(i=1:1:n_concs)
    errorbarxy(exp_mt_lengths, exp_endtags(i,:), exp_errs_x, exp_errs_y(i, :), err_config(i, :));
    line(exp_line_x, exp_line_y(i, :), 'Color', colors(i, :));
end

%title('k on xlink = 0.05 (nM*s)^{-1}; c eff teth = 150 nM; kD teth = 0.2 nM');
xlabel('Length of microtubule (microns)', 'FontSize', 14);
ylabel('Endtag length (microns)', 'FontSize', 14);
%XTick('FontSize', 14);
ylim([0 7]);
xlim([0 16]);