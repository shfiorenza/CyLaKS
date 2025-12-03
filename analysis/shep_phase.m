
clear variables;

site_SID = 0;
xlink_SID = 1;
motor_SID = 2;
endtag_boundary = 0; %500; %125;
flux_win_size = 500;
vel_limit_upper = 10000;
vel_limit_lower = 0;

xlink_analysis_only = false;


seeds = [0];
%{
name = 'motorLife';
sim_name_base = 'shep_0.1nM_%gnM_8_1000_0.6kT_3x_5x_0_motor_%gx';
file_dir = '../out_final_motorLifetime';
output_folder = 'plots_motorLifetime_0_noVelLimit';
var1 = [10, 30, 100];
var1Label = 'Motor concentration (nM)';
var2 = [0.1, 0.3, 1, 3, 10];
var2Label = 'Relative motor lifetime';
%}
%{
name = 'motorVel';
%sim_name_base = 'shep_0.0nM_0.1nM_%i_5000_0.6kT_0_motorMotility_%gx_%gx';
sim_name_base = 'shep_0.1nM_10nM_8_1000_0.6kT_3x_5x_0_motorVelWeighted_%gx_%gx';
file_dir = '../out_final_motorVelWeighted2';
%file_dir = '../out_final_motorMotility';
output_folder = 'plots_motorVelocityWeighted2_0_noVelLimit';
%output_folder = 'plots_motorMotility';
var1 = [0.1, 0.3, 1, 3, 10]; %, 30];
var1Label = 'Relative hydrolysis rate';
var2 = [0.1, 0.3, 1, 3, 10, 30];
var2Label = 'Relative motor lifetime';
%n_pfs = [8, 8, 8, 3, 3, 3];
%}

%{
name = 'xlinkLife';
sim_name_base = 'shep_%gnM_100nM_8_1000_0.6kT_3x_5x_0_xlink_%gx';
file_dir = '../out_final_xlinkLifetime';
output_folder = 'yerp'; %'plots_xlinkLifetime_0_noVelLimit';
var1 = [0.01, 0.03, 0.1, 0.3, 1];
var1Label = 'Crosslinker concentration (nM)';
var2 = [0.1, 0.3, 1, 3, 10];
var2Label = 'Relative crosslinker lifetime';
%}

name = 'xlinkDiffNorm';%'xlinkDiff';
%sim_name_base = 'shep_0.1nM_10nM_8_1000_0.6kT_3x_5x_0_xlinkDiffNorm_%gx_%gx';
sim_name_base = 'shep_0.1nM_10nM_8_1000_0.6kT_0_xlinkDiffNorm_%gx_%gx';
file_dir = '../out_final_xlinkDiffusionNorm';
output_folder = 'plots_xlinkDiffusionNorm_0_noVelLimit';
var1 = [0.0, 0.03, 0.1, 0.3, 1, 3, 10, 30];
%var1 = [0.0, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30];
var1Label = 'Relative lateral diffusion';
var2 = [0.03, 0.1, 0.3, 1, 3, 10, 30];
var2Label = 'Relative longitudinal diffusion';
%}


name = 'protoNum';
sim_name_base = 'shep_0.1nM_50nM_%i_5000_%.1fkT_3x_5x_0';
file_dir = '../out_final_protoNum_5x';
output_folder = 'plots_protoNumber_5x_50nM_0_noVelLimit';
var1 = [1];
var1Label = "";
var2 = [1, 2, 3, 5, 8];
var2Label = "Protofilament Number";
energies = [1.2, 0.8, 0.6, 0.6, 0.6];
%}

run_avg = zeros(length(var1), length(var2), 2);
run_sem = zeros(length(var1), length(var2), 2);
time_avg = zeros(length(var1), length(var2), 2);
time_sem = zeros(length(var1), length(var2), 2);
vel_avg = zeros(length(var1), length(var2), 2);
vel_sem = zeros(length(var1), length(var2), 2);
vel_geo_avg = zeros(length(var1), length(var2), 2);
vel_geo_sem = zeros(length(var1), length(var2), 2);
flux_plus = zeros(length(var1), length(var2), 2);
flux_minus = zeros(length(var1), length(var2), 2);
n_bound_avg = zeros(length(var1), length(var2), 2);

%{
run_theory = zeros(length(var1), length(var2));
vel_theory = zeros(length(var1), length(var2));
time_theory = zeros(length(var1), length(var2));
for i_var = 1 : length(var1)
    for j_var = 1 : length(var2)
            run_theory(i_var, j_var) = 179 * 0.0082 * var2(j_var);
            vel_theory(i_var, j_var) = 8.2 * (1.0/(0.000520 + (0.01333/var1(i_var))));
            time_theory(i_var, j_var) = 179 * var2(j_var) * (0.0009004 + (0.01333/var1(i_var)));
    end
end
%}

for i_var = 1 : length(var1)
    for j_var = 1 : length(var2)
        % for motor + xlink lifetimes + motor vel
        %sim_name = sprintf(sim_name_base, var1(i_var), var2(j_var)) %, seeds(i_seed));
        %sim_name = sprintf(sim_name_base, n_pfs(j_var), var1(i_var), var2(j_var)) %, seeds(i_seed));
        % for xlink diffusion
         % if i_var == 1
         %     sim_name = sprintf("shep_0.1nM_10nM_8_1000_0.6kT_0_xlinkDiffNorm_%gx_0.0x", var2(j_var))
         % else
         %     sim_name = sprintf(sim_name_base, var2(j_var), var1(i_var)) %, seeds(i_seed));
         % end
        % for protoNumer
        sim_name = sprintf(sim_name_base, var2(j_var), energies(j_var))
        try
            params = load_parameters(sprintf('%s/%s', file_dir, sim_name));
            catch_triggered = false;
        catch
            params = load_parameters(sprintf('%s/%s', file_dir, sprintf(sim_name_base, 5, energies(j_var))));
            catch_triggered = true;
        end
        if catch_triggered
            params.n_mts = 8;
            disp("Catch triggered!");
        end
        n_seeds = length(seeds);
 
        % protein ID is unique; make following arrays serial w/ ID as index
        active_proteins = zeros([2 params.n_mts params.n_mts * params.max_sites * 10]);
        n_active = zeros([2 params.n_mts]);
        starting_site = zeros([2 100 * params.n_mts * params.max_sites]) - 1;
        starting_mt = zeros([2 100 * params.n_mts * params.max_sites]) - 1;
        starting_datapoint = zeros([2 100 * params.n_mts * params.max_sites]) - 1;
        runlengths = zeros([2 100 * n_seeds * params.n_mts * params.max_sites]);
        lifetimes = zeros([2 100 * n_seeds  * params.n_mts * params.max_sites]);
        velocities = zeros([2 100 * n_seeds  * params.n_mts * params.max_sites]);
        n_runs = zeros([2 1]);

        midpoint = params.max_sites/2.0;
        for i_seed = 1:n_seeds
            if n_seeds > 1
                disp("FIX SEEDS FIRST!!")
                return
            end
            %sim_name = sprintf('%s_%i', sim_name_base, seeds(i_seed))
            proteinFileStruct = '%s_protein_id.file';
            proteinFileName = sprintf("%s/%s", file_dir, sprintf(proteinFileStruct, sim_name));
            protein_data_file = fopen(proteinFileName);
            if catch_triggered
                raw_motor_data = fread(protein_data_file, '*int');
                size_output = params.max_sites*params.n_mts;
                len_data = size(raw_motor_data);
                n_whole_outputs = floor(len_data/size_output);
                n_whole_outputs = n_whole_outputs(1);
                raw_motor_data = raw_motor_data(1:size_output*n_whole_outputs);
                protein_data = reshape(raw_motor_data, params.max_sites, params.n_mts, n_whole_outputs);
            else
                raw_motor_data = fread(protein_data_file, params.n_mts * params.max_sites * params.n_datapoints, '*int');
                protein_data = reshape(raw_motor_data, params.max_sites, params.n_mts, params.n_datapoints);
            end
            fclose(protein_data_file);

            occuFileStruct = '%s_occupancy.file';
            occuFileName = sprintf("%s/%s", file_dir, sprintf(occuFileStruct, sim_name));
            occu_data_file = fopen(occuFileName);
            if catch_triggered
                raw_occu_data = fread(occu_data_file, '*int');
                raw_occu_data = raw_occu_data(1:size_output*n_whole_outputs);
                occupancy_data = reshape(raw_occu_data, params.max_sites, params.n_mts, n_whole_outputs);
                params.n_datapoints = n_whole_outputs;
            else
                raw_occu_data = fread(occu_data_file, params.n_mts * params.max_sites * params.n_datapoints, '*int');
                occupancy_data = reshape(raw_occu_data, params.max_sites, params.n_mts, params.n_datapoints);
            end
            fclose(occu_data_file);    

            equil_factor = 0.5;
            for i_data = 1 : params.n_datapoints - 1
                if i_data >= params.n_datapoints*equil_factor
                    for i_species = 1 : 1 : 2
                        protein_IDs = protein_data(:, :, i_data);
                        protein_IDs(occupancy_data(:, :, i_data) ~= i_species) = 0;
                        protein_IDs(occupancy_data(:, :, i_data) == i_species) = 1;
                        n_bound = sum(sum(protein_IDs));
                        n_bound_avg(i_var, j_var, i_species) = n_bound_avg(i_var, j_var, i_species) + n_bound/((1.0-equil_factor)*params.n_datapoints);
                    end
                end
                for i_mt = 1: 1 : params.n_mts
                    protein_IDs = protein_data(:, i_mt, i_data);
                    % Scan through IDs of bound motors (-1 means no motor on that site)
                    for i_site = 1 : 1 : params.max_sites
                        i_species = occupancy_data(i_site, i_mt, i_data);
                        if i_species == site_SID || (xlink_analysis_only && i_species == motor_SID)
                            continue;
                        end
                        protein_ID = protein_IDs(i_site);
                        if protein_ID <= 0
                            disp("err");
                            return
                        end
                        % Always count protein on first site
                        if i_site == params.max_sites
                            % Record the protein's starting site if this is the first time
                            % seeing it (-1 means it was not seen last datapoint)
                            if starting_site(i_species, protein_ID) == -1
                                starting_site(i_species, protein_ID) = i_site;
                                starting_datapoint(i_species, protein_ID) = i_data;
                                starting_mt(i_species, protein_ID) =  i_mt;
                                n_active(i_species, i_mt) = n_active(i_species, i_mt) + 1;
                                active_proteins(i_species, i_mt, n_active(i_species, i_mt)) = protein_ID;
                            end
                        % Otherwise if a motor is found, only count first
                        elseif protein_IDs(i_site + 1) ~= protein_ID
                            % calculate flux
                            if i_site >= midpoint - flux_win_size && i_site <= midpoint + flux_win_size
                                for j_mt = i_mt : 1 : params.n_mts + i_mt - 1
                                    if j_mt > params.n_mts
                                        j_mt_adj = j_mt - params.n_mts;
                                    else
                                        j_mt_adj = j_mt;
                                    end
                                    future_site = find(protein_data(:, j_mt_adj, i_data + 1 ) == protein_ID, 1);
                                    if ~isempty(future_site)
                                        if (i_site < midpoint && i_site > midpoint - flux_win_size) && (future_site < midpoint + flux_win_size && future_site > midpoint)
                                            flux_minus(i_var, j_var, i_species) = flux_minus(i_var, j_var, i_species) + 1.0/(params.t_run);
                                        end
                                        if (i_site < midpoint + flux_win_size && i_site > midpoint) && (future_site < midpoint && future_site > midpoint - flux_win_size)
                                            flux_plus(i_var, j_var, i_species) = flux_plus(i_var, j_var, i_species) + 1.0/(params.t_run);
                                        end
                                        break;
                                    end
                                end
                            end
                            % Record the motor's starting site if this is the first time
                            % seeing it (-1 means it was not seen last datapoint)
                            if starting_site(i_species, protein_ID) == -1
                                starting_site(i_species, protein_ID) = i_site;
                                starting_datapoint(i_species, protein_ID) = i_data;
                                starting_mt(i_species, protein_ID) =  i_mt;
                                n_active(i_species, i_mt) = n_active(i_species, i_mt) + 1;
                                active_proteins(i_species, i_mt, n_active(i_species, i_mt)) = protein_ID;
                            end
                        end
                    end
                    % Check one datapoint into the future to see if any motors unbound
                    for i_species = 1 : 1 : 2
                        n_deleted = 0;
                        for i_protein = 1:1:n_active(i_species, i_mt)
                            i_adj = i_protein - n_deleted;
                            protein_ID = active_proteins(i_species, i_mt, i_adj);
                            unbound = true;
                            for j_mt = i_mt : 1 : params.n_mts + i_mt - 1
                                if j_mt > params.n_mts
                                    j_mt_adj = j_mt - params.n_mts;
                                else
                                    j_mt_adj = j_mt;
                                end
                                future_site = find(protein_data(:, j_mt_adj, i_data + 1 ) == protein_ID, 1);
                                if ~isempty(future_site)
                                    % motors have high turnover rate; make sure they do not unbind+rebind in single step
                                    if i_species == motor_SID
                                        if j_mt == starting_mt(i_species, protein_ID) && future_site(1) < starting_site(i_species, protein_ID)
                                            unbound = false;
                                        end
                                    else
                                        unbound = false;
                                    end
                                end
                            end
                            if unbound
                                % Calculate distance traveled
                                end_site = -1;
                                end_mt = -1;
                                for j_mt = i_mt : 1 : params.n_mts + i_mt - 1
                                    if j_mt > params.n_mts
                                        j_mt_adj = j_mt - params.n_mts;
                                    else
                                        j_mt_adj = j_mt;
                                    end
                                    end_site = find(protein_data(:, j_mt_adj, i_data) == protein_ID, 1);
                                    if ~isempty(end_site)
                                        end_mt = j_mt_adj;
                                        break;
                                    end
                                end
                                if isempty(end_site)
                                    disp('Error in finding end site');
                                    return;
                                end
                                start_site = starting_site(i_species, protein_ID);
                                delta = end_site(1) - start_site; % we want actual displacement, not MSD
                                run_length = delta * params.site_size;
                                % Calculate time bound
                                start_datapoint = starting_datapoint(i_species, protein_ID);
                                delta_time = i_data - start_datapoint;
                                run_time = delta_time * params.time_per_datapoint;
                                velocity = (run_length / run_time) * 1000; % convert to nm/s
                                % If time bound is above time cutoff, add to data
                                if end_site(1) > endtag_boundary && run_time >= params.time_per_datapoint && abs(velocity) > vel_limit_lower && abs(velocity) < vel_limit_upper
                                    n_runs(i_species) = n_runs(i_species) + 1;
                                    runlengths(i_species, n_runs(i_species)) = run_length;
                                    lifetimes(i_species, n_runs(i_species)) = run_time;
                                    velocities(i_species, n_runs(i_species)) = velocity;
                                end
                                starting_site(i_species, protein_ID) = -1;
                                starting_datapoint(i_species, protein_ID) = -1;
                                starting_mt(i_species, protein_ID) = -1;
                                % Switch now-deleted entry with last entry in active_motors
                                active_proteins(i_species, i_mt, i_adj) = active_proteins(i_species, i_mt, n_active(i_species, i_mt));
                                active_proteins(i_species, i_mt, n_active(i_species, i_mt)) = -1;
                                n_active(i_species, i_mt) = n_active(i_species, i_mt) - 1;
                                n_deleted = n_deleted + 1;
                            end
                        end
                    end
                end
            end
        end
        % get statistics and generate histograms
        [run_avg(i_var, j_var, xlink_SID), run_sem(i_var, j_var, xlink_SID)] = ...
            get_stats(xlink_SID, runlengths, n_runs, 'normal', "Run length (um)", output_folder, sim_name, "run");
        [time_avg(i_var, j_var, xlink_SID), time_sem(i_var, j_var, xlink_SID)] = ...
            get_stats(xlink_SID, lifetimes, n_runs, 'exponential', "Lifetime (s)", output_folder, sim_name, "time");
        [vel_avg(i_var, j_var, xlink_SID), vel_sem(i_var, j_var, xlink_SID)] = ...
            get_stats(xlink_SID, velocities, n_runs, 'normal', "Velocity (nm/s)", output_folder, sim_name, "vel");
        [vel_geo_avg(i_var, j_var, xlink_SID), vel_geo_sem(i_var, j_var, xlink_SID)] = ...
            get_stats(xlink_SID, abs(velocities), n_runs, 'lognormal', "Geometric velocity (nm/s)", output_folder, sim_name, "velGeo");
        [run_avg(i_var, j_var, motor_SID), run_sem(i_var, j_var, motor_SID)] = ...
            get_stats(motor_SID, abs(runlengths), n_runs, 'exponential', "Run length (um)", output_folder, sim_name, "run");
        [time_avg(i_var, j_var, motor_SID), time_sem(i_var, j_var, motor_SID)] = ...
            get_stats(motor_SID, lifetimes, n_runs, 'exponential', "Lifetime (s)", output_folder, sim_name, "time");
        [vel_avg(i_var, j_var, motor_SID), vel_sem(i_var, j_var, motor_SID)] = ...
            get_stats(motor_SID, velocities, n_runs, 'normal', "Velocity (nm/s)", output_folder, sim_name, "vel");
        %fprintf("(Xlinks) Run: %.2f ± %.2f | Time: %.2f ± %.2f | Vel: %.2f ± %.2f\n", avg_run_xlinks, err_run_xlinks, avg_time_xlinks, err_time_xlinks, avg_vel_xlinks, err_vel_xlinks)
        %fprintf("(Motors) Run: %.2f ± %.2f | Time: %.2f ± %.2f | Vel: %.2f ± %.2f\n", avg_run_motors, err_run_motors, avg_time_motors, err_time_motors, avg_vel_motors, err_vel_motors)

    end
end
%}
net_flux = flux_plus - flux_minus;
if strcmp(name,'protoNum') == 1
    disp("Normalizing flux values for protofilament number!");
    for i = 1 : length(var2)
       net_flux(1, i, :) = net_flux(1, i, :) * 8 / var2(i);
       n_bound_avg(1, i, :) = n_bound_avg(1, i, :) * 8 / var2(i);
    end
end

% phase_runTh = figure('Position', [50, 100, 720, 540]);
% hm = heatmap(run_theory, 'ColorLimits', [0 max(run_theory, [],'all')], 'FontName', 'Arial');
% hm.YDisplayData=flip(hm.YDisplayData);
% phase_velTh = figure('Position', [50, 100, 720, 540]);
% hm = heatmap(vel_theory, 'ColorLimits', [0 max(vel_theory, [],'all')], 'FontName', 'Arial');
% hm.YDisplayData=flip(hm.YDisplayData);
% phase_timeTh = figure('Position', [50, 100, 720, 540]);
% hm = heatmap(time_theory, 'ColorLimits', [0 max(time_theory, [],'all')], 'FontName', 'Arial');
% hm.YDisplayData=flip(hm.YDisplayData);

for SID =  1 : 1 : 2
    make_phase(abs(run_avg), var1, var2, var1Label, var2Label, "displacement (um)", output_folder, "run", name, SID)
    make_phase(time_avg, var1, var2, var1Label, var2Label, "lifetime (s)", output_folder, "time", name, SID)
    make_phase(abs(vel_avg), var1, var2, var1Label, var2Label, "velocity (nm/s)", output_folder, "vel", name, SID)
    make_phase(net_flux, var1, var2, var1Label, var2Label, "flux (1/s)", output_folder, "flux", name, SID)
    make_phase(n_bound_avg, var1, var2, var1Label, var2Label, "number bound", output_folder, "num", name, SID)
    make_plot(abs(run_avg), var1, var2, var1Label, var2Label, "displacement (um)", output_folder, "run", name, SID);
    make_plot(time_avg, var1, var2, var1Label, var2Label, "lifetime (s)", output_folder, "time", name, SID);
    make_plot(abs(vel_avg), var1, var2, var1Label, var2Label, "velocity (nm/s)", output_folder, "vel", name, SID);
    make_plot(net_flux, var1, var2, var1Label, var2Label, "flux (1/s)", output_folder, "flux", name, SID);
    make_plot(n_bound_avg, var1, var2, var1Label, var2Label, "number bound", output_folder, "num", name, SID);
    if SID == xlink_SID
        make_phase(abs(vel_geo_avg), var1, var2, var1Label, var2Label, "geometric velocity (nm/s)", output_folder, "velGeo", name, SID)
        make_plot(abs(vel_geo_avg), var1, var2, var1Label, var2Label, "geometric velocity (nm/s)", output_folder, "velGeo", name, SID);
    end
end

function [avg, sem] = get_stats(SID, data_arr, n_runs_arr, distribution, dataLabel, output_folder, sim_name, plotname)
n_runs = n_runs_arr(SID);
data = data_arr(SID, 1:n_runs);
if strcmp(distribution, 'lognormal') == 1
    avg = exp(sum(log(data)) / n_runs);
    sd = exp(std(log(data)));
    sem = sd / sqrt(n_runs);
else
    avg = sum(data) / n_runs;
    sd = std(data);
    sem = sd / sqrt(n_runs);
end
if n_runs < 10
    disp("Not enough statistics for histogram fitting")
    return
end
if strcmp(distribution, 'exponential') == 1
    min_val = min(data);
    data = data - min_val; 
end
pd = fitdist(data', distribution);
conf_inv = paramci(pd);
if strcmp(distribution, 'normal') == 1
    mu = pd.mu;
    sigma = pd.sigma;
    stderr = (conf_inv(2) - conf_inv(1)) / (2*1.96); % get std_err
elseif strcmp(distribution, 'exponential') == 1
    mu = pd.mu + min_val;
    sigma = nan;
    stderr = (conf_inv(2) - conf_inv(1)) / (2*1.96); % get std_err
elseif strcmp(distribution, 'lognormal') == 1
    mu = exp(pd.mu);
    sigma = exp(pd.sigma);
    stderr = (exp(conf_inv(2)) - exp(conf_inv(1))) / (2*1.96); % get std_err
else
    disp("unrecognized probability distribution")
    return
end
if abs(avg - mu) > 1e-6
    fprintf("Error calculating %s for species %i\n", dataLabel, SID)
    avg = nan;
    sem = nan;
    return;
end
switch plotname
    case 'run'
        bin_lower_lim = 0.1;
    case 'time'
        bin_lower_lim = 0.2;
    case 'vel'
        bin_lower_lim = 50;
    case 'velGeo'
        bin_lower_lim = 50;
    otherwise 
        disp('error')
        return
end
if SID == 1
    speciesLabel = "xlink";
elseif SID == 2
    speciesLabel = "motor";
end
fig = figure('Position', [50 50 720 600]);
set(gcf, 'DefaultAxesFontName', 'Arial');
set(gcf, 'DefaultTextFontName', 'Arial');
bin_width = max(bin_lower_lim, 2*iqr(data)/(n_runs^(1/3)));
n_bins = ceil((max(data) - min(data))/bin_width);
h = histfit(data, n_bins, distribution);
h(1).FaceAlpha = 0.5;
h(1).EdgeAlpha = 0.5;
h(1).FaceColor = [0 0 0];
h(1).EdgeColor = 'none';
h(1).BarWidth = 0.75;
h(2).LineWidth = 3;
xlim([min(0, min(data)) max(data)])
xlabel(sprintf("%s of %ss", dataLabel, speciesLabel));
ylabel("Count");
dim1 = [0.625 0.7 0.2 0.2];
str1 = sprintf('Avg = %#.3g ± %#.2g', avg, sem);
dim2 = [0.625 0.65 0.2 0.2];
str2 = sprintf('Std = %#.3g', sd);
dim3 = [0.625 0.6 0.2 0.2];
str3 = sprintf('Mu = %#.3g ± %#.2g', mu, stderr);
dim4 = [0.625 0.55 0.2 0.2];
str4 = sprintf('Sigma = %#.3g', sigma);
dim5 = [0.625 0.5 0.2 0.2];
str5 = sprintf('N = %i', n_runs);
dim6 = [0.625 0.45 0.2 0.2];
str6 = sprintf('Bin width = %#.3g', bin_width);
annotation('textbox', dim1, 'String', str1, 'FitBoxToText', 'on', 'EdgeColor','none', 'FontSize', 14, "FontName", "Arial");
annotation('textbox', dim2, 'String', str2, 'FitBoxToText', 'on', 'EdgeColor','none', 'FontSize', 14, "FontName", "Arial");
annotation('textbox', dim3, 'String', str3, 'FitBoxToText', 'on', 'EdgeColor','none', 'FontSize', 14, "FontName", "Arial");
annotation('textbox', dim4, 'String', str4, 'FitBoxToText', 'on', 'EdgeColor','none', 'FontSize', 14, "FontName", "Arial");
annotation('textbox', dim5, 'String', str5, 'FitBoxToText', 'on', 'EdgeColor','none', 'FontSize', 14, "FontName", "Arial");
annotation('textbox', dim6, 'String', str6, 'FitBoxToText', 'on', 'EdgeColor','none', 'FontSize', 14, "FontName", "Arial");
set(gca,'box','off')
set(gca, 'FontSize', 18);
set(gca, 'FontName', 'Arial');
set(gca,'TickDir','out');
set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);
title(sprintf("%s (%ss)", sim_name, speciesLabel),'Interpreter', 'none', "FontName", "Arial", "FontSize", 10);
exportgraphics(fig, sprintf('%s/%s_%s_%s.png', output_folder, sim_name, speciesLabel, plotname),'ContentType','image');
exportgraphics(fig, sprintf('%s/%s_%s_%s.pdf', output_folder, sim_name, speciesLabel, plotname),'ContentType','vector');
%saveas(fig, sprintf('%s/%s_%s_%s.svg', output_folder, sim_name, speciesLabel, plotname),'svg');
%saveas(fig, sprintf('%s/%s_%s_%s.png', output_folder, sim_name, speciesLabel, plotname),'png');
close(fig);
%fprintf("%.3g +/- %.1g\n", mu, stderr);
%fprintf("%.3g +/- %.1g\n", avg, sem);
end

function make_phase(data, var1, var2, var1Label, var2Label, dataLabel, output_folder, plotname, name, SID)
data = data(:, :, SID);
if SID == 1
    speciesLabel = "MAP";
elseif SID == 2
    speciesLabel = "motor";
else
    disp("Error saving files.")
    return;
end
fig = figure('Position', [50, 100, 720, 540]);
hm = heatmap(data, 'ColorLimits', [0 max(data, [],'all')], 'FontName', 'Arial');
hm.YDisplayLabels = num2str( var1' );
hm.YLabel = var1Label; 
hm.XDisplayLabels = num2str( var2' );
hm.XLabel = var2Label;
hm.CellLabelFormat = '%.3g';
hm.YDisplayData=flip(hm.YDisplayData);
hm.GridVisible = 'off';
hm.Units = 'centimeters'; 
hm.Position(3:4) = min(hm.Position(3:4))*[1,1]; % make heatmap square
set(gca, "FontSize", 16);
set(gca, "FontName", "Arial");
ylabel(hm.NodeChildren(2), "Average " + speciesLabel + " " + dataLabel, "FontSize", 18, "FontName", "Arial");
saveas(fig, sprintf('%s/phase_%s_%s_%s.png', output_folder, speciesLabel, plotname, name), 'png');
saveas(fig, sprintf('%s/phase_%s_%s_%s.svg', output_folder, speciesLabel, plotname, name), 'svg');
end

function make_plot(data, var1, var2, var1Label, var2Label, dataLabel, output_folder, plotname, name, SID)
data = data(:, :, SID);
if SID == 1
    speciesLabel = "MAP";
elseif SID == 2
    speciesLabel = "motor";
else
    disp("Error saving files.")
    return;
end
fig = figure('Position', [50 50 960 540]);
if strcmp(name,'protoNum') == 1
    plot([1 2 3 5 8], data', '.', 'MarkerSize', 50)
    xticks([1 2 3 5 8])
    xlim([0 10]);
else
    plot(data', '.', 'MarkerSize', 50)
    xticks(1:length(var2));
    xticklabels(num2str( var2' ));
    xlim([0 length(var2)+1]);
end
ylim([0 1.15*max(data)]);
xlabel(var2Label);
ylabel("Average " + speciesLabel + " " + dataLabel);
legendLabel = num2str( var1' );
leg = legend(legendLabel, 'Location', 'northeastoutside');
title(leg, var1Label);
set(gca, 'FontName', 'Arial')
set(gca, 'FontSize', 18);
set(gca,'box','off')
set(gca,'TickDir','out');
set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);
pbaspect([1 1 1]);
saveas(fig, sprintf('%s/plot_%s_%s_%s.png', output_folder, speciesLabel, plotname, name), 'png');
saveas(fig, sprintf('%s/plot_%s_%s_%s.svg', output_folder, speciesLabel, plotname, name), 'svg');
end

