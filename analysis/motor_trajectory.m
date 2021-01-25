clear all
n_mts = 1;
n_sites = 2000;
n_datapoints = 10000;
t_per_datapoint = 0.1;

ID = 1435; % 462; %622;

fileDirectory = '/home/shane/projects/CyLaKS/%s';
fileName = 'hybrid_motor_0.05_0_protein_id.file';
data_file = fopen(sprintf(fileDirectory, fileName));
mt_array = fread(data_file, [n_mts*n_sites, n_datapoints], '*int');
fclose(data_file);

trajectory = zeros([n_datapoints 1]) - 10;
first_sighting = true;
for i_data = 1 : n_datapoints
    for i_site = 1 : n_sites
        if mt_array(i_site, i_data) == ID
            trajectory(i_data, 1) = i_site;
            if first_sighting
                n_start = i_data;
                first_sighting = false;
            end
            n_end = i_data;
            break;
        end
    end
end

%%plot fig%%
fig1 = figure(1);
set(fig1,'Position', [50, 50, 2*480, 2*270])

% reverse t so as to mirror a real kymograph
x = trajectory(n_start : n_end);
t =  linspace(n_end, n_start, n_end - n_start + 1);
plot(x, t, ':ob', 'LineWidth', 1.5, 'MarkerSize', 5);

xlim([0 n_sites]);
ylim([n_start-1 n_end+1]);
xlabel('Site on microtubule');
ylabel({'Datapoint number', sprintf('[Motor ID: %d]', ID)});
%yticks(linspace(n_start, n_end, 5));
yticklabels({'1000', '900', '800', '700', '600','500','400','300','200','100'});