clear all;
sim_name = 'xlink_bind_ii_test';
kbT = 4.114;                    % pN * nm
k_on = 0.000238;                % 1/(nM*s)
c_bind = 4500;                  % nM
k_off_ii = 14.3;                % 1/s
dist_cutoff = 4;                % no. sites
mt_coords = [0, 0];
n_datapoints = 1000000;
% 'reported' probabilities are from sim's built-in probability tracker
p_bind_reported = [2.69e-05, 2.61e-05, 1.81e-05, 4.61e-06, 2.16e-07];
p_unbind_reported = [0.000357, 0.000369, 0.000529, 0.00208, 0.0450];

speciesID = 1;
static_head_pos = dist_cutoff;
mt_lengths = [2*dist_cutoff + 1, 2*dist_cutoff + 1];
n_mts = length(mt_lengths);
max_length = max(mt_lengths);

fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
fileStruct = '%s_occupancy.file';
data_file = fopen(sprintf(fileDirectory, sprintf(fileStruct, sim_name)));
raw_data = fread(data_file, max_length * n_mts * n_datapoints, '*int');
occu_data = reshape(raw_data, max_length, n_mts, n_datapoints);
fclose(data_file);

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

% Generate figure
fig1 = figure();
set(fig1,'Position', [50, 50, 2.5*480, 2.5*300])
% Plot avg_occu data
plot(linspace(-dist_cutoff, dist_cutoff, 2*dist_cutoff + 1), ... 
    avg_occu, '*', 'MarkerSize', 14, 'LineWidth', 3);
hold on
% Plot theoretical avg_occu curve for continuous space using raw params
%   (note: uses delta_u function defined at the end of this script)
k_d = k_off_ii / k_on;

avg_occu_theory = @(x) c_bind / (c_bind + k_d*exp(delta_u(x)/kbT));
fplot(@(x) avg_occu_theory(x), 'LineWidth', 2);

%{
% Plot theoretical avg_occu at each point using p values reported by sim
for x_dist = 0 : dist_cutoff
    p_bind = p_bind_reported(x_dist + 1);
    p_unbind = p_unbind_reported(x_dist + 1);
    avg_occu_reported = p_bind / (p_bind + p_unbind);
    plot(x_dist, avg_occu_reported, 'or', 'MarkerSize', 20, 'LineWidth', 2);
    plot(-x_dist, avg_occu_reported, 'or', 'MarkerSize', 20, 'LineWidth', 2);
end
%}
% Label them axes
xlabel('Crosslinker extension (n\_sites)');
ylabel('Average occupancy');
% Generate legend
legend(["Simulation data", "Analytic theory"]);
% Force avg_occu plot on top of avg_occu_theory plot
h = get(gca,'Children');
set(gca,'Children',[h(2) h(1)])

% delta_u function used to get spring energy for avg_occu_theory
function u = delta_u(x)
    k_spring = 0.453; % pN / nm
    site_size = 8; % nm
    dr = sqrt(32^2 + (x*site_size)^2) - 32;
    u = 0.5*k_spring*dr*dr;
end
