
clear variables;

%sim_name_base = 'shep_%inM_%inM_%i_%skT';
%file_dir = '../shepherding_baseline'; 

%sim_name_base = 'shep_%gnM_%gnM_8_1000_0.6kT_3x_5x_%i';
%file_dir = '../out_final';

%sim_name_base = 'shep_0.1nM_%gnM_8_1000_0.6kT_3x_5x_0_motor_%gx';
%file_dir = '../out_final_motorLifetime';

%sim_name_base = 'shep_%gnM_100nM_8_1000_0.6kT_3x_5x_0_xlink_%gx';
%file_dir = '../out_final_xlinkLifetime';

sim_name_base = 'shep_0.1nM_10nM_8_1000_0.6kT_3x_5x_0_xlinkDiff_%gx_%gx';
file_dir = '../out_final_xlinkDiffusion4';

%sim_name_base = 'shep_0.1nM_50nM_%i_1000_%.1fkT_3x_5x_0';
%file_dir = '../outputProto5';

%const1 = 10;
%const1Label = 'nM PRC1';
%const2 = 200;
%const2Label = 'nM K401';

%var1 = [0.1, 1];
%var1Label = 'Crosslinker concentration (nM)';
%var2 = [1, 10, 50, 100, 250, 500, 1000];
%var2Label = 'Motor concentration (nM)';
%seeds = [0, 1, 2, 3, 4, 5];

%var1 = [10, 30, 100];
%var1Label = 'Motor concentration (nM)';
%var2 = [0.1, 0.3, 1, 3, 10];
%var2Label = 'Relative motor lifetime';

%var1 = [0.03, 0.1, 0.3, 1]; %[0.01, 0.03, 0.1, 0.3, 1];
%var1Label = 'Crosslinker concentration (nM)';
%var2 = [0.1, 0.3, 1, 3, 10];
%var2Label = 'Relative crosslinker lifetime';
seeds = [0];

var1 = [0.1, 0.3];%, 1, 3, 10];
var1Label = 'Relative longitudinal diffusion';
var2 = [0.1, 0.3];%, 1, 3, 10];
var2Label = 'Relative lateral diffusion';

%var1 = [1];
%var1Label = "";
%var2 = [1, 2, 3, 5, 8];
%var2Label = "Protofilament Number";
%energies = [1.2, 0.8, 0.6, 0.6, 0.6];
%seeds = [0];

titleOccu = "Average Total Occupancy";%sprintf('Occupancy for %i %s and %i %s', const1, const1Label, const2, const2Label);
titleRatio = "Average Occupancy Ratio";%sprintf('Plus- vs minus-end occupancy for %i %s and %i %s', const1, const1Label, const2, const2Label);

start_frame = 3000;

ratio_window = 50; % in n_sites

% number of crosslinkers bound to MT at each timestep
occu_avg = zeros(length(var1), length(var2));
occu_err = zeros(length(var1), length(var2));

occu_plus = zeros(length(var1), length(var2), length(seeds));
occu_minus = zeros(length(var1), length(var2), length(seeds));

% ratio of plus-end to minus-end occupancy at each timestep
ratio_avg = zeros(length(var1), length(var2));
ratio_err = zeros(length(var1), length(var2));

site_ID = 0;
xlink_ID = 1;
motor_ID = 2;

for i_var = 1 : length(var1)
    for j_var = 1 : length(var2)
        %if i_var == 1 && j_var > 4
        %    continue
        %end
        for i_seed = 1 : length(seeds)
            %sim_name = sprintf(sim_name_base, const1, const2, var1(i_var), var2(j_var));
            sim_name = sprintf(sim_name_base, var1(i_var), var2(j_var)); %, seeds(i_seed));
            %sim_name = sprintf(sim_name_base, var2(j_var), energies(j_var));

            % Load parameter structure
            params = load_parameters(sprintf('%s/%s', file_dir, sim_name));
            % Open occupancy data file
            occupancy_filename = sprintf('%s/%s_occupancy.file', file_dir, sim_name);
            occupancy = zeros(params.max_sites, params.n_mts, params.n_datapoints) - 1;
            occupancy = load_data(occupancy, occupancy_filename, '*int');

            n_sites_tot = 0;
            for i_mt = 1 : params.n_mts
                n_sites_tot = n_sites_tot + params.mt_lengths(i_mt);
            end
            n_active_frames = params.n_datapoints - start_frame;
            for i_data = start_frame : params.n_datapoints
                for i_mt = 1 : params.n_mts
                    for i_site = 1 : params.mt_lengths(i_mt)
                        sid = occupancy(i_site, i_mt, i_data);  % species label ID
                        if sid == xlink_ID
                            occu_avg(i_var,j_var) = occu_avg(i_var,j_var) + 1/(n_active_frames * n_sites_tot);
                            if i_site <= ratio_window
                                occu_plus(i_var, j_var, i_seed) = occu_plus(i_var, j_var, i_seed) + 1/(n_active_frames * n_sites_tot);
                            elseif i_site > params.mt_lengths(i_mt) - ratio_window
                                occu_minus(i_var, j_var, i_seed) = occu_minus(i_var, j_var, i_seed) + 1/(n_active_frames * n_sites_tot);
                            end
                        elseif sid == motor_ID
                            % nothing yet
                        end
                    end
                end

                %{
            ratio = occu_plus / occu_minus;
            if isnan(ratio)
                disp(occu_plus)
                disp(occu_minus)
                disp('bro');
            end
            ratio_avg(i_var, j_var) = ratio_avg(i_var, j_var) + ratio/n_active_frames;
                %}
            end
        end
    end
end
%}

%{
for i_var = 1 : length(var1)
    for j_var = 1 : length(var2)
        ratio_avg(i_var, j_var) = occu_plus(i_var, j_var) / occu_minus(i_var, j_var);
    end
end
%}

ratio_avg = mean(occu_plus, 3) ./ mean(occu_minus, 3);

fig1 = figure('Position', [50, 100, 540, 540]);
hm = heatmap(occu_avg);
hm.Title = titleOccu;
hm.YDisplayLabels = string(var1);
hm.YLabel = var1Label; 
hm.XDisplayLabels = string(var2);
hm.XLabel = var2Label;
hm.CellLabelFormat = '%.3f';
hm.YDisplayData=flip(hm.YDisplayData);
%clim([0 1])
set(gca, 'FontSize', 18);

fig2 = figure('Position', [590, 100, 540, 540]);
hm = heatmap(ratio_avg);
hm.Title = titleRatio;
hm.YDisplayLabels = string(var1);
hm.YLabel = var1Label; 
hm.XDisplayLabels = string(var2);
hm.XLabel = var2Label;
hm.CellLabelFormat = '%.1f';
hm.YDisplayData=flip(hm.YDisplayData);
clim([0 100])
set(gca, 'FontSize', 18);

%fig3 = figure('Position', [50 50 720 600]);
fig3 = figure('Position', [50 50 1000 600]);
%semilogy(var2, ratio_avg, '.', 'MarkerSize', 50)
%semilogx(var2, ratio_avg, '.', 'MarkerSize',50)
%loglog(var2, ratio_avg, '.', 'MarkerSize',50)
semilogx(var1, ratio_avg, '.', 'MarkerSize',50)

set(gca,'box','off')
set(gca, 'FontSize', 28);
set(gca,'TickDir','out');
set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);

ylabel("Tip occupancy ratio");
%yticks([0 250 500])
%xlabel(var2Label)
xlabel(var1Label)
%xlim([0 10]);
xlim([0 12]);
%xticks(var2);
xticks(var1)
%legendLabel = {'0.03 nM', '0.1 nM', '0.3 nM', '1.0 nM'};
%legendLabel = {'10 nM', '30 nM', '100 nM'};
legendLabel = {'0.1', '0.3', '1', '3', '10'};
legend(legendLabel, 'Location', 'northeastoutside');