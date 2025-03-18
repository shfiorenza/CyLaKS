clear variables;

%sim_name = 'out_endtags1/shep_0.1nM_100nM_8_500_0.2kT_1x_0';
sim_name = 'output16/shep_1nM_100nM_8_0.2kT_1x_0';
sim_name = 'output18/shep_1nM_100nM_8_1000_0.8kT_1x_0';
sim_name = 'shep_50x_0.02_0.5kT_0.131_0.131_0.1nM_10nM';
sim_name = 'shep_10x_0.01_1kT_0.131_0.131_1nM_100nM';
sim_name = 'output22/shep_1nM_100nM_8_250_0.2kT_0.1x_0.3x_1';
sim_name = 'output26/shep_1nM_100nM_8_1000_0.6kT_1.75x_0.3x_0';
sim_name = 'output23/shep_1nM_100nM_8_1000_0.6kT_0.1x_1.2x_0';
sim_name = 'output24/shep_1nM_100nM_8_1000_0.6kT_2x_1.5x_0';
sim_name = 'output28/shep_1nM_100nM_8_1000_0.6kT_3x_5x_0';
sim_name = 'out_final/shep_0.1nM_100nM_8_1000_0.6kT_3x_5x_0';
sim_name = 'out_final_newCombos/shep_0.75nM_30nM_8_1000_0.6kT_3x_5x_0';
 

% Load parameter structure
file_dir = '..'; 
params = load_parameters(sprintf('%s/%s', file_dir, sim_name));

% Open occupancy data file 
occupancy_filename = sprintf('%s/%s_occupancy.file', file_dir, sim_name);
occupancy = zeros(params.max_sites, params.n_mts, params.n_datapoints) - 1;
occupancy = load_data(occupancy, occupancy_filename, '*int');

xlink_speciesID = 1;
motor_speciesID = 2;

xlink_raw_data = occupancy; 
motor_raw_data = occupancy; 

xlink_raw_data(xlink_raw_data ~= xlink_speciesID) = 0;
xlink_raw_data(xlink_raw_data == xlink_speciesID) = 1;
motor_raw_data(motor_raw_data ~= motor_speciesID) = 0;
motor_raw_data(motor_raw_data == motor_speciesID) = 1;

xlink_avg_occupancy = zeros([params.n_datapoints params.n_mts]);
motor_avg_occupancy = zeros([params.n_datapoints params.n_mts]);

motor_avg_occupancy_tot = zeros([params.n_datapoints 1]);
xlink_avg_occupancy_tot = zeros([params.n_datapoints 1]);

for i_data = 1 : params.n_datapoints
    for i_mt = 1 : params.n_mts
        n_sites = params.mt_lengths(i_mt);
        for i_site = 1 : n_sites
            species_id = occupancy(i_site, i_mt, i_data);
            if species_id == 0
                continue;
            elseif species_id == xlink_speciesID
                xlink_avg_occupancy(i_data, i_mt) = xlink_avg_occupancy(i_data, i_mt) + 1.0 / double(n_sites);
                xlink_avg_occupancy_tot(i_data, 1) = xlink_avg_occupancy_tot(i_data, 1) + 1.0 / double(n_sites * params.n_mts);
            elseif species_id == motor_speciesID
                motor_avg_occupancy(i_data, i_mt) = motor_avg_occupancy(i_data, i_mt) + 1.0 / double(n_sites);
                motor_avg_occupancy_tot(i_data, 1) = motor_avg_occupancy_tot(i_data, 1) + 1.0 / double(n_sites * params.n_mts);
            end
        end
    end
end

color_xlink = [12 220 210] / 255;
color_motor = [214 77 156] / 255;

fig_xlink = figure('Position',[50 50 720 540]);
plot(linspace(0, 0.1 * params.n_datapoints, params.n_datapoints), xlink_avg_occupancy_tot, ... 
    'LineWidth', 5, 'Color', color_xlink);
xlabel("Time (s)");
%xlim([0 200]);
ylabel("Fractional occupancy");
%ylim([0 0.05]);
set(gca,'box','off')
set(gca, 'FontSize', 24);

fig_motor = figure('Position',[50 50 720 540]);
plot(linspace(0, 0.1 * params.n_datapoints, params.n_datapoints), motor_avg_occupancy_tot, ...
    'LineWidth', 5, 'Color', color_motor);
xlabel("Time (s)");
%xlim([0 200]);
ylabel("Fractional occupancy");
%ylim([0.15 0.2]);
set(gca,'box','off')
set(gca, 'FontSize', 24);

saveas(fig_xlink, 'occu_vs_t_xlink.png', 'png');
saveas(fig_motor, 'occu_vs_t_motor.png', 'png');

