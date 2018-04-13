clear all
n_datapoints = 100000;
length_of_microtubule = 250;
species_ID = 2;
k_off_frac = 1;
k_on = 0.015;

temp_one = zeros([length_of_microtubule 1]);
final_mt = zeros([length_of_microtubule 1]);

polarityArray = {'Plus-end on right'};

%fileDirectory = '/media/shane/Shane''s External HDD (1 TB)/Parameter Scan 1/%s';
fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
fileName = 'test21_occupancy.file';
%fileName = '5_50.000_100.00_1250_occupancy.file'

data_file = fopen(sprintf(fileDirectory, fileName));
raw_data = fread(data_file, [length_of_microtubule, n_datapoints], '*int');
fclose(data_file);

raw_data((raw_data ~= species_ID) | (raw_data == 0)) = 0;
raw_data((raw_data == species_ID) | (raw_data ~= 0)) = 1;

% Read in and average occupancy data over all datapoints
for i=1:1:n_datapoints
    temp_one(:, 1) = raw_data(:, i);
    final_mt(:, 1) = final_mt(:, 1) + double(temp_one(:, 1)./n_datapoints);
end


smoothed_final_mt = smoothdata(final_mt);
smoothed_final_mt = final_mt;


% Read through MT to find approx end-tag position (where the intensity
% [really occupancy in our case] is 0.5 of it's maxmium [1])
% alt method: highest slope
highest_slope = 0;
endtag_site = length_of_microtubule;
for i=1:1:length_of_microtubule
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

    if(site_occupancy > 0.5)
        endtag_site =  i;
        break;
    end
end

%%plot fig%%
fig1 = figure(1);
set(fig1,'Position', [50, 50, 2.5*480, 2.5*300])
plot(linspace(0, length_of_microtubule*0.008, length_of_microtubule), smoothed_final_mt);
% Put vertical red line where endtag starting position is
hold on
endtag_pos = endtag_site * 0.008;
plot([endtag_pos endtag_pos], [0 1], ':r', 'LineWidth', 0.1);

%%style stuff%%
grid on
grid minor
title(sprintf('%g micron endtag for %g micron MT\n k on = %g s^-^1, k off (stalled) = k off / %i', ...
    (length_of_microtubule - endtag_site) * 0.008, length_of_microtubule * 0.008, k_on, k_off_frac));

xlabel({'Distance along microtubule (microns)'}); %, sprintf('(%d microns)', ...
    %length_of_microtubule * 8 / 1000)});
ylabel('Fraction of the time occupied');
axis = gca;
axis.YLim = [0 1];
axis.XLim = [0 length_of_microtubule*0.008];
axis.TickDir = 'out';
axis.Box = 'off';
axis.GridLineStyle = '-';
set(findall(axis, 'Type', 'Line'), 'LineWidth', 2);
legend(polarityArray, 'Location', 'northeastoutside');
