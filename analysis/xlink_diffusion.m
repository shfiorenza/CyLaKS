
clear variables;
fileDirectory = '/home/shane/projects/CyLaKS/%s';
sim_name_base = 'test_xlink_bind_ii';
sim_name_base = 'test_xlink_diffusion_Boltzmann';
offsets = [0]; %[0.0, 4.1]; % , 8.2];
labels = [""]; % ["0.0", "4.1"]; %, "8.2"];
dist_cutoff = 8;                % no. sites

normalize_nonRest = 0;
k_spring = 0.0453;  % pN/ nm

kbT = 4.114;                    % pN * nm
k_on = 0.000238;                % 1/(nM*s)
c_bind = 5000;                  % nM
k_off_ii = 0.143;                % 1/s
speciesID = 1;
k_d = k_off_ii / k_on;
static_head_pos = dist_cutoff;
mt_lengths = [2*dist_cutoff + 1, 2*dist_cutoff + 1];
n_mts = length(mt_lengths);
max_length = max(mt_lengths);
n_offsets = length(offsets);

avg_occu = zeros(n_offsets, 2 * dist_cutoff + 1);
for i_offset = 1 : n_offsets
    % Occupancy data -- species ID of each occupant or -1 for none
    if n_offsets > 1
        sim_name = sprintf('%s_%s', sim_name_base, labels(i_offset));
    else 
        sim_name = sim_name_base;
    end
    log_file = sprintf(fileDirectory, sprintf('%s.log', sim_name));
    log = textscan(fileread(log_file), '%s %s', 'Delimiter', '=');
    params = log{1, 1};
    values = log{1, 2};
    n_datapoints = str2double(values{contains(params, "n_datapoints ")});
    
    occu_file = fopen(sprintf(fileDirectory, sprintf('%s_occupancy.file', sim_name)));
    raw_occu_data = fread(occu_file, max_length * n_mts * n_datapoints, '*int');
    fclose(occu_file);
    occu_data = reshape(raw_occu_data, max_length, n_mts, n_datapoints);
    
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
            disp('this should not happen')
        else
            % Get x_dist for this configuration (-dist_cutoff to +dist_cutoff)
            x_dist = dynamic_head_pos - static_head_pos;
            % Convert to an index between 1 and 2*dist_cutoff + 1
            index = x_dist + dist_cutoff;
            % All x_dists but this one should not factor this step into avg_occu
            n_points_applicable(index) = n_points_applicable(index) + 1;
        end
    end
    % Calculate the average occupancy for each x_dist
    %for i_data = 1 : n_datapoints
     %   dynamic_head_pos = find(occu_data(:, 2, i_data));
        %if ~isempty(dynamic_head_pos)
            % Get x_dist for this configuration (-dist_cutoff to +dist_cutoff)
            %x_dist = dynamic_head_pos - static_head_pos;
            % Convert to an index between 1 and 2*dist_cutoff + 1
            %index = x_dist + dist_cutoff;
            % All x_dists but this one should not factor this step into avg_occu
            for index = 1 : 2*dist_cutoff+1
                if normalize_nonRest
                    avg_occu(i_offset, index) = n_points_applicable(index) / n_points_applicable(dist_cutoff);
                else
                    avg_occu(i_offset, index) = n_points_applicable(index) / n_points_applicable(dist_cutoff+1);
                end
            end
            %avg_occu(i_offset, index) = avg_occu(i_offset, index) + double(1) / n_points_applicable(index);
        %end
    %end
end
%}
% Generate figure
fig1 = figure();
set(fig1,'Position', [50, 50, 720, 540])
hold all;
% Plot avg_occu data
color = [0 0.447 0.741; 0.85, 0.325, 0.098; 0.929, 0.694, 0.125; ...
    0.494, 0.184, 0.556; 0.466, 0.674, 0.188; 0.301, 0.745, 0.933];
for i_offset = 1 : n_offsets
        % Plot theoretical avg_occu curve for continuous space using raw params
    %   (note: uses delta_u function defined at the end of this script)
    avg_occu_theory = @(x) exp(-delta_u(x, offsets(i_offset), k_spring)/kbT);
    %avg_occu_theory = @(x) c_bind / (c_bind + k_d*exp(delta_u(x, offsets(i_offset))/kbT));
    fplot(@(x) avg_occu_theory(x), 'LineWidth', 2, 'Color', [0.6 0.6 0.6]);
    sim_data = plot(linspace(-dist_cutoff, dist_cutoff, 2*dist_cutoff + 1), ...
        avg_occu(i_offset, :), 'o', 'MarkerSize', 12,'MarkerEdgeColor', color(i_offset, :), ...
        'LineWidth', 3);
    sim_data.MarkerFaceColor = sim_data.MarkerEdgeColor;
    sim_data.Color = sim_data.MarkerFaceColor;

end

%xlim([-(dist_cutoff+1) dist_cutoff+1]);
%xticks([-4 -2 0 2 4]);
%ylim([-0.025 0.75]);
%yticks([0 0.2 0.4 0.6]);



% Label them axes
xlabel('Displacement along lattice', 'FontSize', 22);
ylabel('Average occupancy', 'FontSize', 22);
set(gca, 'FontSize', 22);
% Generate legend

%fracLabels = ["0", "2", "4", "6"];
%legendLabel = fracLabels + "-nm offset";
%{
legendLabel = ["Simulation (0-nm offset)", "Simulation (4-nm offset)", "Theory"];
h = get(gca,'Children');
legend([h(3) h(1) h(4)], legendLabel, 'FontSize', 20, 'location', ... 
    'northwest', 'NumColumns', 1);
legend('boxoff');
%}

%new_handle = copyobj(hlegend1,);                 %copy legend 1 --> retain
%legend(hplot2, 'Data 2', 'Location','West');            %display legend 2
% Force avg_occu plot on top of avg_occu_theory plot

%set(gca,'Children',[h(2) h(1)])

function u = delta_u(x, offset, k_spring)
    site_size = 8.2; % nm
    dr = sqrt(32^2 + (x*site_size + offset)^2) - 32;
    u = 0.5*k_spring*dr*dr;
end