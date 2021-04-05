
clear variables;
sim_name_base = 'mobility_both';
concs = [20,50,80,120,220,420];
seeds = [0,1,2,3];

n_mts = 1;
n_sites = 1250;
site_size = 0.008;

n_concs = length(concs);
n_seeds = length(seeds);
fileDirectory = '../run_mobility_both/%s';

counts = zeros(n_concs, n_seeds, n_sites);
for i_conc = 1 : n_concs
    for i_seed = 1 : n_seeds
        sim_name = sprintf('%s_%i_%i', sim_name_base, concs(i_conc), seeds(i_seed))
        % Open log file and parse it into param labels & their values
        log_file = sprintf(fileDirectory, sprintf('%s.log', sim_name));
        log = textscan(fileread(log_file), '%s %s', 'Delimiter', '=');
        params = log{1, 1};
        values = log{1, 2};
        % Read in system params
        n_datapoints = str2double(values{contains(params, "n_datapoints ")});
        % Use actual recorded number of datapoints to parse thru data/etc
        if any(contains(params, "N_DATAPOINTS ") ~= 0)
            n_datapoints = str2double(values{contains(params, "N_DATAPOINTS ")});
        end
        proteinFileName = '%s_motorID.file';
        proteinFile = sprintf(fileDirectory, sprintf(proteinFileName, sim_name));
        file = fopen(proteinFile);
        data = fread(file, n_mts * n_sites * n_datapoints, '*int');
        fclose(file);
        protein_ids = reshape(data, n_sites, n_mts, n_datapoints);
        for i_data = 1 : n_datapoints
            for i_site = 1 : n_sites
                id = protein_ids(i_site, 1, i_data);
                if id ~= -1
                    for j_site = i_site : n_sites
                        neighb_id = protein_ids(j_site, 1, i_data);
                        if neighb_id ~= -1 && neighb_id ~= id
                            dist = j_site - i_site;
                            counts(i_conc, i_seed, dist) = counts(i_conc, i_seed, dist) + 1;
                            i_site = j_site; 
                        end
                    end
                end
            end
        end
    end
end
%}

avg_counts = zeros(n_concs, n_sites);
err_counts = zeros(n_concs, n_sites);

for i_conc = 1: n_concs
    for i_site = 1 : n_sites
       avg_counts(i_conc, i_site) = mean(counts(i_conc, :, i_site)); 
    end
   avg_counts(i_conc, :) = avg_counts(i_conc, :) / max(avg_counts(i_conc, 2:n_sites));
end

fig1 = figure;
set(fig1, 'Position', [50 50 1000 600]);
hold all;
smoothed_counts = zeros(n_concs, n_sites - 1);
for i_conc = 1 : n_concs
    smoothed_counts(i_conc, :) = smoothdata(avg_counts(i_conc, 2:n_sites), 'gaussian', 30);
    plot(linspace(1 * site_size, n_sites * site_size, n_sites - 1), ... 
        smoothed_counts(i_conc, :), 'LineWidth', 2)
end
ylim([0 1]);
ylabel("Normalized counts");
xlabel("Distance between motors (microns)");
legend(concs + " pM", 'Location', 'northeast', 'FontSize', 22)
set(gca, 'FontSize', 24);

avg_distance = zeros(n_concs, 1);
for i_conc = 1 : n_concs
    %pd = fitdist(avg_counts(i_conc, 2:n_sites)','exponential');
    pd = fitdist(avg_counts(i_conc, :)','exponential')
    avg_distance(i_conc) = pd.mu;
end

fig2 = figure();
set(fig2, 'Position', [75 75 1000 600]);
hold all;
plot(concs, avg_distance);

fig3 = figure();
h = histfit(smoothed_counts(1, 2:1249), 100, 'exponential');