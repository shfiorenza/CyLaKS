function params = get_parameters(sim_name)

% Open log file and parse it into param labels & their values
fileDirectory = '../%s';
%fileDirectory = '/home/shane/data_kif4a_paper/run_mobility_both/%s';
log_file = sprintf(fileDirectory, sprintf('%s.log', sim_name));
log = textscan(fileread(log_file), '%s %s', 'Delimiter', '=');
params = log{1, 1};
values = log{1, 2};
% Read in system params
dt = sscanf(values{contains(params, "dt ")}, '%g');
time_per_datapoint = sscanf(values{contains(params, "t_snapshot ")}, '%g');
n_datapoints = str2double(values{contains(params, "n_datapoints ")});
% Use actual recorded number of datapoints to parse thru data/etc
if any(contains(params, "N_DATAPOINTS ") ~= 0)
    n_datapoints = str2double(values{contains(params, "N_DATAPOINTS ")});
end
site_size =  sscanf(values{contains(params, "site_size ")}, '%g') / 1000; % in um
% Read in number of MTs
n_mts = sscanf(values{contains(params, "count ")}, '%g');
if any(contains(params, "COUNT ") ~= 0)
    n_mts = sscanf(values{contains(params, "COUNT ")}, '%g');
end
if any(contains(params, "n_subfilaments") ~= 0)
    n_sub = sscanf(values{contains(params, "n_subfilaments ")}, '%g');
    if n_sub > n_mts
       n_mts = n_sub;
    end
end
% Read in MT lengths (in n_sites)
mt_lengths = zeros(1, n_mts);
for i_mt = 1 : n_mts
    string = sprintf("n_sites[%i] ", i_mt - 1);
    mt_lengths(i_mt) =  sscanf(values{contains(params, string)}, '%i');
    if any(contains(params, sprintf("N_SITES[%i] ", i_mt - 1)) ~= 0)
        string = sprintf("N_SITES[%i] ", i_mt - 1);
        mt_lengths(i_mt) = sscanf(values{contains(params, string)}, '%i');
    end
end
max_sites = max(mt_lengths);
polarity = zeros(1, n_mts);
for i_mt = 1 : n_mts
    string = sprintf("polarity[%i] ", i_mt - 1);
    polarity(i_mt) =  sscanf(values{contains(params, string)}, '%i');
end


n_dims = 2; % hard-coded for now; CyLaKS always outputs data in 2-D
xlink_cutoff = 5; % FIXME: dynamically get this from log
teth_cutoff = 10; % FIXME: dynamically get this from log


end