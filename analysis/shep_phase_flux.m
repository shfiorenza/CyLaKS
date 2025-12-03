% FIX ME -- very janky implementation of multi PFs at the moment

clear variables;

seeds = [0];
%{
name = 'motorLife';
sim_name_base = 'shep_0.1nM_%gnM_8_1000_0.6kT_3x_5x_0_motor_%gx';
file_dir = '../out_final_motorLifetime';
output_folder = 'plots_flux';
var1 = [10, 30, 100];
var1Label = 'Motor concentration (nM)';
var2 = [0.1, 0.3, 1, 3, 10];
var2Label = 'Relative motor lifetime';
%}

name = 'motorVel';
sim_name_base = 'shep_0.1nM_10nM_8_1000_0.6kT_0_motorVelConstNum_%gx_%gx';
%sim_name_base = 'shep_0.1nM_10nM_8_1000_0.6kT_3x_5x_0_motorVel_%gx_%gx';
file_dir = '../out_final_motorVel_constNum';
output_folder = 'plots_motorVelocityConstNum';
%file_dir = '../out_final_motor_velocity';
%output_folder = 'plots_motorVelocity';
var1 = [0.1, 0.3, 1, 3, 10, 30];
var1Label = 'Relative hydrolysis rate';
var2 = [0.1, 0.3, 1, 3, 10, 30];
var2Label = 'Relative motor lifetime';
%}

%{
name = 'xlinkLife';
sim_name_base = 'shep_%gnM_100nM_8_1000_0.6kT_3x_5x_0_xlink_%gx';
file_dir = '../out_final_xlinkLifetime';
output_folder = 'plots_flux';
var1 = [0.03, 0.1, 0.3, 1];
var1Label = 'Crosslinker concentration (nM)';
var2 = [0.1, 0.3, 1, 3, 10];
var2Label = 'Relative crosslinker lifetime';
%}
%{
name = 'xlinkDiff';
%sim_name_base = 'shep_0.1nM_10nM_8_1000_0.6kT_3x_5x_0_xlinkDiff_%gx_%gx';
sim_name_base = 'shep_0.1nM_10nM_8_1000_0.6kT_0_xlinkDiffNorm_%gx_%gx';
file_dir = '../out_final_xlinkDiffusionNorm';
output_folder = 'plots_xlinkDiffusionNorm';
var1 = [0.0, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30];
var1Label = 'Relative lateral diffusion';
var2 = [0.03, 0.1, 0.3, 1, 3, 10, 30];
var2Label = 'Relative longitudinal diffusion';
%}
%{
name = 'protoNum';
sim_name_base = 'shep_0.1nM_50nM_%i_1000_%.1fkT_3x_5x_0';
file_dir = '../out_final_proto';
output_folder = 'plots_flux';
var1 = [1];
var1Label = "";
var2 = [1, 2, 3, 5, 8];
var2Label = "Protofilament Number";
energies = [1.2, 0.8, 0.6, 0.6, 0.6];
%}


xlink_SID = 1;
chosen_SID = xlink_SID;

flux_plus = zeros(length(var1), length(var2));
flux_minus = zeros(length(var1), length(var2));

for i_var = 1 : length(var1)
    for j_var = 1 : length(var2)
        % for motor + xlink lifetimes + motor vel
        sim_name = sprintf(sim_name_base, var1(i_var), var2(j_var)) %, seeds(i_seed));
        % for xlink diffusion
        % if i_var == 1
        %     sim_name = sprintf("shep_0.1nM_10nM_8_1000_0.6kT_0_xlinkDiffNorm_%gx_0.0x", var2(j_var))
        % else
        %     sim_name = sprintf(sim_name_base, var2(j_var), var1(i_var)) %, seeds(i_seed));
        % end
        % for protoNumer
        %sim_name = sprintf(sim_name_base, var2(j_var), energies(j_var))

        % Open log file and parse it into param labels & their values
        %log_file = sprintf('%s/%s', file_dir, sprintf('%s_0.log', sim_name_base));
        log_file = sprintf('%s/%s', file_dir, sprintf('%s.log', sim_name));
        log = textscan(fileread(log_file), '%s %s', 'Delimiter', '=');
        params = log{1, 1};
        values = log{1, 2};
        % Read in number of MTs
        n_mts = sscanf(values{contains(params, "count ")}, '%g');
        if any(contains(params, "n_subfilaments") ~= 0)
            n_sub = sscanf(values{contains(params, "n_subfilaments ")}, '%g');
            if n_sub > n_mts
                n_mts = n_sub;
            end
        end
        if any(contains(params, "COUNT ") ~= 0)
            n_mts = sscanf(values{contains(params, "COUNT ")}, '%g');
        end
        mt_lengths = zeros(1, n_mts);
        for i_mt = 1 : n_mts
            string = sprintf("n_sites[%i] ", i_mt - 1);
            mt_lengths(i_mt) = sscanf(values{contains(params, string)}, '%i');
            if any(contains(params, sprintf("N_SITES[%i] ", i_mt - 1)) ~= 0)
                string = sprintf("N_SITES[%i] ", i_mt - 1);
                mt_lengths(i_mt) = sscanf(values{contains(params, string)}, '%i');
            end
        end
        n_sites = max(mt_lengths);
        % Read in system params
        t_tot = sscanf(values{contains(params, "t_run ")}, '%g');
        dt = sscanf(values{contains(params, "dt ")}, '%g');
        steps_per_datapoint = str2double(values{contains(params, "n_steps_per_snapshot ")});
        time_per_datapoint = dt * steps_per_datapoint;
        n_datapoints = str2double(values{contains(params, "n_datapoints ")});
        % Use actual recorded number of datapoints to parse thru data/etc
        if any(contains(params, "N_DATAPOINTS ") ~= 0)
            n_datapoints = str2double(values{contains(params, "N_DATAPOINTS ")});
        end
        n_dims = 2;
        site_size = 0.0082; % in um
        n_seeds = length(seeds);



        for i_seed = 1:n_seeds
            if n_seeds > 1
                disp("FIX SEEDS FIRST!!")
                return
            end
            %sim_name = sprintf('%s_%i', sim_name_base, seeds(i_seed))
            proteinFileStruct = '%s_protein_id.file';
            proteinFileName = sprintf("%s/%s", file_dir, sprintf(proteinFileStruct, sim_name));
            protein_data_file = fopen(proteinFileName);
            raw_motor_data = fread(protein_data_file, n_mts * n_sites * n_datapoints, '*int');
            fclose(protein_data_file);
            protein_data = reshape(raw_motor_data, n_sites, n_mts, n_datapoints);
            occuFileStruct = '%s_occupancy.file';
            occuFileName = sprintf("%s/%s", file_dir, sprintf(occuFileStruct, sim_name));
            occu_data_file = fopen(occuFileName);
            raw_occu_data = fread(occu_data_file, n_mts * n_sites * n_datapoints, '*int');
            fclose(occu_data_file);
            occupancy_data = reshape(raw_occu_data, n_sites, n_mts, n_datapoints);
            for i_data = 1 : n_datapoints - 1
                for i_mt = 1 : 1 : n_mts
                    protein_IDs = protein_data(:, i_mt, i_data);
                    % Scan through IDs of bound motors (-1 means no motor on that site)
                    for i_site = midpoint - win_size:1:midpoint + win_size
                        occupant = occupancy_data(i_site, i_mt, i_data);
                        if occupant ~= chosen_SID
                            continue;
                        end
                        protein_ID = protein_IDs(i_site);
                        % if looking at motors, only check one head
                        if protein_ID > 0 && protein_IDs(i_site - 1) ~= protein_ID
                            start_site = i_site;
                            for j_mt = 1 : 1 : n_mts
                                future_site = find(protein_data(:, j_mt, i_data + 1 ) == protein_ID, 1);
                                if ~isempty(future_site)
                                    if start_site < midpoint && future_site > midpoint
                                        flux_minus(i_var, j_var) = flux_minus(i_var, j_var) + 1.0/(t_tot);
                                    end
                                    if start_site > midpoint && future_site < midpoint
                                        flux_plus(i_var, j_var) = flux_plus(i_var, j_var) + 1.0/(t_tot);
                                    end
                                    break;
                                end
                            end
                        end
                    end
                end
            end
        end
        disp(flux_plus(i_var, j_var));
        disp(flux_minus(i_var, j_var));
    end
end
net_flux = flux_plus - flux_minus
%}
if strcmp(name,'protoNum') == 1
    for i = 1 : length(var2)
       net_flux(i) = net_flux(i) * 8 / var2(i);
    end
end

phase_flux = figure('Position', [50, 100, 720, 540]);
hm = heatmap(net_flux, 'ColorLimits', [0 max(net_flux, [],'all')], 'FontName', 'Arial');
hm.YDisplayLabels = num2str( var1' );
hm.YLabel = var1Label; 
hm.XDisplayLabels = num2str( var2' );
hm.XLabel = var2Label;
hm.CellLabelFormat = '%.2f';
hm.YDisplayData=flip(hm.YDisplayData);
hm.GridVisible = 'off';
hs = struct(hm);
set(gca, 'FontSize', 16);
ylabel(hs.Colorbar, "Average MAP flux (1/s)", 'FontSize', 20);



plot_flux = figure('Position', [50 50 720 540]);
%plot_flux = figure('Position', [50 50 720 600]);
set(gcf, 'DefaultAxesFontName', 'Arial');
set(gcf, 'DefaultTextFontName', 'Arial');
plot(net_flux', '.', 'MarkerSize', 50)
%plot([1 2 3 5 8],net_flux', '.', 'MarkerSize', 50)
xticks(1:length(var2));
xticklabels(num2str( var2' ));
xlabel(var2Label);
ylim([0 inf]);
ylabel("Average MAP flux (1/s)");
legendLabel = num2str( var1' );
legend(legendLabel, 'Location', 'northeastoutside');
set(gca,'box','off')
set(gca, 'FontSize', 18);
%set(gca, 'FontSize', 28);
set(gca,'TickDir','out');
set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);
xlim([0 10]);
%xticks([1 2 3 5 8])
%}

saveas(phase_flux, sprintf('%s/phase_flux_%s.svg', output_folder, name), 'svg');
saveas(phase_flux, sprintf('%s/phase_flux_%s.png', output_folder, name), 'png');
saveas(plot_flux, sprintf('%s/plot_flux_%s.svg', output_folder, name), 'svg');
saveas(plot_flux, sprintf('%s/plot_flux_%s.png', output_folder, name), 'png');
%close all