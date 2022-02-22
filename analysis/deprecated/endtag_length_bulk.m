clear variables;
xlink_concs = [0]; % , 1, 4];
base_name = "endtagE_%i_%i_%i";
mt_lengths = [2, 4, 6, 8, 10, 14]; % in microns
seeds = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
% experimental parameters
exp_mt_lengths = [2.4, 4.0, 5.6, 7.2, 8.8, 10.3, 13.5];
exp_endtags_0 = [1.2, 1.3, 1.4, 1.6, 1.65, 1.75, 2.25];
exp_errs_0_y = [0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.2];
exp_endtags_1 = [1.4, 1.6, 1.7, 2.0, 2.1, 2.4, -1];
exp_errs_1_y = [0.1, 0.6, 0.2, 0.3, 0.2, 0.3, 0];
exp_endtags_4 = [1.6, 2.1, 2.4, 3.1, 3.6, 4.0, 5.3];
exp_errs_4_y = [0.2, 0.6, 0.25, 0.4, 0.25, 1.0, 0.25];
exp_endtags = [exp_endtags_0.', exp_endtags_1.', exp_endtags_4.'].';
exp_errs_x = [0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8];
exp_errs_y = [exp_errs_0_y.', exp_errs_1_y.', exp_errs_4_y.'].';
exp_line_x = [2 14.5];
exp_line_y = [[1.15 2.23]; [1.37 2.88]; [1.67 5.05]];
motor_speciesID = 2;
xlink_speciesID = 1;
n_datapoints = 10000;
starting_point = 5000;
active_datapoints = n_datapoints - starting_point;
fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
fileStruct = '%s_occupancy.file';

n_concs = length(xlink_concs);
n_lengths = length(mt_lengths);
n_seeds = length(seeds);
endtag_data = zeros(n_concs, n_lengths, n_seeds);

for i_conc = 1:1:n_concs
    xlink_conc = xlink_concs(i_conc);

    for i_length = 1:1:n_lengths

        for i_seed = 1:1:n_seeds
            n_sites = mt_lengths(i_length) * 125;
            simName = sprintf(base_name, xlink_conc, n_sites, seeds(i_seed));
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
            for i = starting_point:1:n_datapoints
                motor_avg_occupancy(:, 1) = motor_avg_occupancy(:, 1) + double(motor_raw_data(:, i)) ./ active_datapoints;
                xlink_avg_occupancy(:, 1) = xlink_avg_occupancy(:, 1) + double(xlink_raw_data(:, i)) ./ active_datapoints;
            end

            smooth_window = n_sites / 50;
            motor_occupancy = smoothdata(motor_avg_occupancy, 'movmean', smooth_window);
            xlink_occupancy = smoothdata(xlink_avg_occupancy, 'movmean', smooth_window);
            net_occupancy = motor_occupancy + xlink_occupancy;
            occupancy_slope = smoothdata(gradient(net_occupancy, 0.008), 'movmean', smooth_window);
            occupancy_accel = smoothdata(gradient(occupancy_slope, 0.008), 'movmean', smooth_window);
            max_occupancy = max(net_occupancy);
            min_slope = min(occupancy_slope);

            past_threshold = false;
            i_threshold = 1;
            endtag_site = 0;

            for i_site = 1:n_sites

                if (~past_threshold && net_occupancy(i_site) < 0.5 * max_occupancy)
                    past_threshold = true;
                    i_threshold = i_site;
                end

                if (past_threshold && occupancy_slope(i_site) > (min_slope / 2))
                    endtag_site = i_site;
                    break;
                end

            end

            if endtag_site == 0
                [max_accel, i_peak] = max(occupancy_accel(i_threshold:n_sites));
                endtag_site = i_threshold + i_peak;
            end

            endtag_length = endtag_site * 0.008;

            endtag_data(i_conc, i_length, i_seed) = endtag_length;
        end

    end

end

avg_endtag_length = zeros(n_concs, n_lengths);
err_endtag_length = zeros(n_concs, n_lengths);

for i_conc = 1:n_concs

    for i_length = 1:n_lengths
        % Get average final length across all seeds
        mean_endtag_length = mean(endtag_data(i_conc, i_length, :));
        avg_endtag_length(i_conc, i_length) = mean_endtag_length;
        % Calculate standard deviation of this average
        variance = 0.0;

        for i_seed = 1:n_seeds
            diff_sq = (mean_endtag_length - endtag_data(i_conc, i_length, i_seed))^2;
            variance = variance + diff_sq / (n_seeds - 1);
        end

        err_endtag_length(i_conc, i_length) = sqrt(variance);
    end

end

%final_data = endtag_data(1,:,:);

set(0, 'DefaultLineLineWidth', 1.25);
%set(0,'DefaultLegendAutoUpdate','off');
fig1 = figure(1);
set(fig1, 'Position', [50, 50, 2 * 480, 2 * 300])
hold on

colors = ['b'; 'r'; 'k'];
err_config = [{'^b', 'b', 'b'}; {'^r', 'r', 'r'}; {'^k', 'k', 'k'}];

for i = 1:1:n_concs
    errorbar(mt_lengths, avg_endtag_length(i_conc, :), err_endtag_length(i_conc, :), ...
        'o', 'LineWidth', 2); % , 'Color', colors(i, :));
    %plot(mt_lengths, endtag_data(1, :, i),'*','LineWidth',2, 'MarkerSize', 15, 'Color', colors(i,:));
end

%legendLabel = cellstr(num2str(xlink_concs' / 10, '%#.1f nM PRC1'));
%legendLabel = {'0.0 nM PRC1', '0.1 nM PRC1', '0.4 nM PRC1'};

for (i = 1:1:n_concs)
    line(exp_line_x, exp_line_y(i, :), 'LineStyle', '--', 'Color', 'red'); %, 'Color', colors(i, :));
end

legend({'Simulation', 'Experiment'});
line([0 16], [0 0], 'Linestyle', ':', 'Color', 'k');

%title('k on xlink = 0.05 (nM*s)^{-1}; c eff teth = 150 nM; kD teth = 0.2 nM');
xlabel('Length of microtubule (microns)', 'FontSize', 14);
ylabel('Endtag length (microns)', 'FontSize', 14);
%XTick('FontSize', 14);
%ylim([0 6]);
%xlim([0 16]);
