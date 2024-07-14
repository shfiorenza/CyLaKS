function params = load_parameters(sim_name)

% Open log file and parse it into param labels & their values
log_file = sprintf('%s.log', sim_name);
log = textscan(fileread(log_file), '%s %s', 'Delimiter', '=');
params = log{1, 1};
values = log{1, 2};
% Read in system params
%wut = values{contains(params, "dt ")}
dt = sscanf(values{contains(params, "dt ")}, '%f');
time_per_datapoint = sscanf(values{contains(params, "t_snapshot ")}, '%g');
n_datapoints = str2double(values{contains(params, "n_datapoints ")});
% Use actual recorded number of datapoints to parse thru data/etc
if any(contains(params, "N_DATAPOINTS ") ~= 0)
    n_datapoints = str2double(values{contains(params, "N_DATAPOINTS ")});
end
site_size =  sscanf(values{contains(params, "site_size ")}, '%g') / 1000; % in um
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

params = struct('dt',dt,'time_per_datapoint', time_per_datapoint, ...
    'n_datapoints', n_datapoints, 'site_size', site_size, 'n_mts', n_mts, ...
    'mt_lengths', mt_lengths, 'max_sites', max_sites, 'polarity', polarity, ...
    'n_dims', n_dims);

end