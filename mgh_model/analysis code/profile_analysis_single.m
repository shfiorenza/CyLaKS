clear all
% Often-changed variables
n_sites = 250;
simName = 'testf'
% Pseudo-constant variables
motor_speciesID = 2;
xlink_speciesID = 1;
n_datapoints = 100000;
starting_point = 50000;
active_datapoints = n_datapoints - starting_point;

temp_one = zeros([n_sites 1]);
final_mt = zeros([n_sites 1]);

polarityArray = {'Plus-end on left'};

%fileDirectory = '/media/shane/Shane''s External HDD (1 TB)/Parameter Scan 1/%s';
%fileDirectory = '/home/shane/Desktop/pseudo_crackpot/%s';
fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
fileStruct = '%s_occupancy.file';
fileName = sprintf(fileDirectory, sprintf(fileStruct, simName));

data_file = fopen(fileName);
raw_data = fread(data_file, [n_sites, n_datapoints], '*int');
fclose(data_file);

raw_data((raw_data ~= motor_speciesID) & (raw_data ~= xlink_speciesID)) = 0;
raw_data((raw_data == motor_speciesID) | (raw_data == xlink_speciesID)) = 1;

% Read in and average occupancy data over all datapoints
for i=starting_point:1:n_datapoints
    temp_one(:, 1) = raw_data(:, i);
    final_mt(:, 1) = final_mt(:, 1) + double(temp_one(:, 1)./active_datapoints);
end


smoothed_final_mt = smoothdata(final_mt);
smoothed_final_mt = final_mt;


% Read through MT to find approx end-tag position (where the intensity
% [really occupancy in our case] is 0.5 of it's maxmium [1])
% alt method: highest slope
highest_slope = 0;
endtag_site = n_sites;
for i=n_sites:-1:1
    site_occupancy = smoothed_final_mt(i, 1);
    if(i == 1)
        prev_site_occupancy = 0;
    else
        prev_site_occupancy = smoothed_final_mt(i - 1, 1);
    end
 %   slope = site_occupancy - prev_site_occupancy;
  %  if(slope > highest_slope && site_occupancy > 0.4)
   %     highest_slope = slope;
    %    endtag_site = i;
    %end 

    if(site_occupancy > 0.7)
        endtag_site =  i;
        break;
    end
end

%%plot fig%%
fig1 = figure(1);
set(fig1,'Position', [50, 50, 2.5*480, 2.5*300])
plot(linspace(0, n_sites*0.008, n_sites), smoothed_final_mt);
% Put vertical red line where endtag starting position is
hold on
endtag_pos = endtag_site * 0.008;
plot([endtag_pos endtag_pos], [0 1], ':r', 'LineWidth', 0.1);

%%style stuff%%
grid on
grid minor
title(sprintf('%g micron endtag for %g micron MT', ... % \n k on = %g s^-^1, k off (stalled) = k off / %i', ...
    (endtag_site) * 0.008, n_sites * 0.008)); % k_on, k_off_frac));

xlabel({'Distance along microtubule (microns)'}); %, sprintf('(%d microns)', ...
    %length_of_microtubule * 8 / 1000)});
ylabel('Fraction of the time occupied');
axis = gca;
axis.YLim = [0 1];
axis.XLim = [0 n_sites*0.008];
axis.TickDir = 'out';
axis.Box = 'off';
axis.GridLineStyle = '-';
set(findall(axis, 'Type', 'Line'), 'LineWidth', 2);
legend(polarityArray, 'Location', 'northeastoutside');
