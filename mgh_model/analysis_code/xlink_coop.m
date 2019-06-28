clear variables;
simName = 'coop1';
mt_length = [1000];
n_mts = length(mt_length);
n_sites_max = max(mt_length);
n_steps = 100000000;
delta_t = 0.00001;
n_datapoints = 10000;
starting_point = 0;
active_datapoints = n_datapoints - starting_point;
xlink_speciesID = 1;

fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
fileStructure = '%s_occupancy.file';

xlink_concs = [1, 4,8];
n_concs = length(xlink_concs);
norm_intensity = zeros([n_datapoints n_concs]);

for i_conc=1:1:n_concs
    simName = sprintf('coop%i', xlink_concs(i_conc));
    
    fileName = sprintf(fileDirectory, sprintf(fileStructure, simName));
    xlink_data_file = fopen(fileName);
    xlink_raw_data = fread(xlink_data_file, [n_mts * n_sites_max * n_datapoints], '*int');
    fclose(xlink_data_file);
    xlink_data = reshape(xlink_raw_data, n_sites_max, n_mts, n_datapoints);
    
    xlink_data(xlink_data ~= xlink_speciesID) = 0;
    xlink_data(xlink_data == xlink_speciesID) = 1;
    
    for i_data=1:1:n_datapoints
        for i_mt=1:1:n_mts
            occupancy = xlink_data(:, i_mt, i_data);
            n_sites = mt_length(i_mt);
            for i_site=1:1:n_sites
                norm_intensity(i_data, i_conc) = norm_intensity(i_data) + double(occupancy(i_site))/n_sites;
            end
        end
    end
end

norm_intensity = smoothdata(norm_intensity);
intensity_avgs = zeros([n_concs 1]);
for i_conc=1:1:n_concs
   intensity_avgs(i_conc) =  mean(norm_intensity(5000:10000, i_conc));
end

plot(xlink_concs, intensity_avgs, '*--', 'LineWidth', 2);