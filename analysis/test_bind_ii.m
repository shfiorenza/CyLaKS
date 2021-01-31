%
clear variables;

sim_name = 'xlink_bind_ii_8.2';
offset = 8.2;

fileDirectory = '/home/shane/projects/CyLaKS/%s';
kbT = 4.114;                    % pN * nm
k_on = 0.000238;                % 1/(nM*s)
c_bind = 500000;                  % nM
k_off_ii = 143;                % 1/s
dist_cutoff = 4;                % no. sites
n_datapoints = 1000000; % 0 % extra factor of 10 for 0.0 case

speciesID = 1;
static_head_pos = dist_cutoff;
mt_lengths = [2*dist_cutoff + 1, 2*dist_cutoff + 1];
n_mts = length(mt_lengths);
max_length = max(mt_lengths);

% Occupancy data -- species ID of each occupant or -1 for none
occu_file = fopen(sprintf(fileDirectory, sprintf('%s_occupancy.file', sim_name)));
raw_occu_data = fread(occu_file, max_length * n_mts * n_datapoints, '*int');
fclose(occu_file);
occu_data = reshape(raw_occu_data, max_length, n_mts, n_datapoints);
%{
% MT position data -- gives coordinate of plus/minus end of MTs
pos_file = fopen(sprintf(fileDirectory, sprintf('%s_filament_pos.file', sim_name)));
raw_pos_data = fread(pos_file, 2*n_dims * n_mts * n_datapoints, '*double');
fclose(pos_file);
filament_pos = reshape(raw_pos_data, n_dims, 2, n_mts, n_datapoints);
%}

occu_data(occu_data ~= speciesID) = 0;
occu_data(occu_data == speciesID) = 1;

% Determine the number of datapoints that each occupancy should be averaged
% over in order to take into account that the head cannot occupy more than
% one site at a time
n_points_applicable = zeros(2 *dist_cutoff + 1, 1);
for i_data = 1 : n_datapoints
   dynamic_head_pos = find(occu_data(:, 2, i_data));
   % If second head isn't bound, datapoint is applicable for all datapoints
   if isempty(dynamic_head_pos)
       n_points_applicable = n_points_applicable + 1;
   else
       % Get x_dist for this configuration (-dist_cutoff to +dist_cutoff)
       x_dist = dynamic_head_pos - static_head_pos;
       % Convert to an index between 1 and 2*dist_cutoff + 1
       index = x_dist + dist_cutoff;
       % All x_dists but this one should not factor this step into avg_occu
       n_points_applicable(index) = n_points_applicable(index) + 1;
   end
end

avg_occu = zeros(2 * dist_cutoff + 1, 1);
% Calculate the average occupancy for each x_dist
for i_data = 1 : n_datapoints
   dynamic_head_pos = find(occu_data(:, 2, i_data));
   if ~isempty(dynamic_head_pos)
       % Get x_dist for this configuration (-dist_cutoff to +dist_cutoff)
       x_dist = dynamic_head_pos - static_head_pos;
       % Convert to an index between 1 and 2*dist_cutoff + 1
       index = x_dist + dist_cutoff;
       % All x_dists but this one should not factor this step into avg_occu
       avg_occu(index) = avg_occu(index) + double(1) / n_points_applicable(index);
   end
end
%}

% Generate figure
fig1 = figure();
set(fig1,'Position', [50, 50, 480, 720])
% Plot avg_occu data
sim_data = plot(linspace(-dist_cutoff, dist_cutoff, 2*dist_cutoff + 1), ... 
    avg_occu, 'o', 'MarkerSize', 14,'MarkerEdgeColor', [0 0.447 0.741], ...
    'LineWidth', 3);
sim_data.MarkerFaceColor = sim_data.MarkerEdgeColor;
sim_data.Color = sim_data.MarkerFaceColor;
hold on
xlim([-dist_cutoff dist_cutoff]);
%xticks
ylim([0 0.65]);
yticks([0 0.2 0.4 0.6]);

% Plot theoretical avg_occu curve for continuous space using raw params
%   (note: uses delta_u function defined at the end of this script)
k_d = k_off_ii / k_on;
avg_occu_theory = @(x) c_bind / (c_bind + k_d*exp(delta_u(x, offset)/kbT));
fplot(@(x) avg_occu_theory(x), 'LineWidth', 2);

% Label them axes
xlabel('Lattice site difference', 'FontSize', 24);
ylabel('Average occupancy', 'FontSize', 24);
set(gca, 'FontSize', 24);
% Generate legend
legend(["Simulation", "Theory"], 'FontSize', 20,...
    'location', 'northwest');
% Force avg_occu plot on top of avg_occu_theory plot
h = get(gca,'Children');
set(gca,'Children',[h(2) h(1)])

% delta_u function used to get spring energy for avg_occu_theory

function u = delta_u(x, offset)
    k_spring = 0.453; % pN / nm
    site_size = 8.2; % nm
    dr = sqrt(32^2 + (x*site_size + offset)^2) - 32;
    u = 0.5*k_spring*dr*dr;
end
