clear variables;

sim_name_base = 'shep_%inM_%inM_%i_%skT';
file_dir = '../shepherding_baseline'; 

const1 = 10;
const1Label = 'nM PRC1';
const2 = 200;
const2Label = 'nM K401';

var1 = [8, 4, 2, 1];
var1Label = 'Number of protofilaments';
var2 = ["0.0", "0.75", "1.5", "2.25"];
var2Label = 'Interaction Energy (kBT)';

titleOccu = sprintf('Occupancy for %i %s and %i %s', const1, const1Label, const2, const2Label);
titleRatio = sprintf('Plus- vs minus-end occupancy for %i %s and %i %s', const1, const1Label, const2, const2Label);

start_frame = 2000;

ratio_window = 100; % in n_sites

% number of crosslinkers bound to MT at each timestep
occu_avg = zeros(length(var1), length(var2));
occu_err = zeros(length(var1), length(var2));

occu_plus = zeros(length(var1), length(var2));
occu_minus = zeros(length(var1), length(var2));

% ratio of plus-end to minus-end occupancy at each timestep
ratio_avg = zeros(length(var1), length(var2));
ratio_err = zeros(length(var1), length(var2));

site_ID = 0;
xlink_ID = 1;
motor_ID = 2;

for i_var = 1 : length(var1)
    for j_var = 1 : length(var2)
        sim_name = sprintf(sim_name_base, const1, const2, var1(i_var), var2(j_var));
        
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
                            occu_plus(i_var, j_var) = occu_plus(i_var, j_var) + 1/(n_active_frames * n_sites_tot);
                        elseif i_site > params.mt_lengths(i_mt) - ratio_window
                            occu_minus(i_var, j_var) = occu_minus(i_var, j_var) + 1/(n_active_frames * n_sites_tot);
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
%}

%{
for i_var = 1 : length(var1)
    for j_var = 1 : length(var2)
        ratio_avg(i_var, j_var) = occu_plus(i_var, j_var) / occu_minus(i_var, j_var);
    end
end
%}

ratio_avg = occu_plus ./ occu_minus;

fig1 = figure();
hm = heatmap(occu_avg);
hm.Title = titleOccu;
hm.YDisplayLabels = string(var1);
hm.YLabel = var1Label; 
hm.XDisplayLabels = string(var2);
hm.XLabel = var2Label;
caxis([0 1]);

fig2 = figure();
hm = heatmap(ratio_avg);
hm.Title = titleRatio;
hm.YDisplayLabels = string(var1);
hm.YLabel = var1Label; 
hm.XDisplayLabels = string(var2);
hm.XLabel = var2Label;