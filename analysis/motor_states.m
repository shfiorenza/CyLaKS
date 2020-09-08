clear all;
% Often-changed variables
n_sites = 50000;
simName = 'processivity_20pM_0';
% Pseudo-constant variables
n_mts = 1;
n_datapoints = 10000;

fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
motorFileStruct = '%s_motorID.file';

fileName = sprintf(fileDirectory, sprintf(motorFileStruct, simName));

data_file = fopen(fileName);
mt_array = fread(data_file, [n_mts*n_sites, n_datapoints], '*int');
fclose(data_file);

n_hits_single = 0;
n_hits_double = 0;

for i_data=1:1:n_datapoints
    mt = mt_array(:, i_data);
    for i_site=1:1:n_sites
        ID = mt(i_site);
        if ID ~= -1
            if i_site == 1
                if mt(i_site + 1) == ID
                    n_hits_double = n_hits_double + 1;
                else
                    n_hits_single = n_hits_single + 1;
                end
            elseif i_site == n_sites
                if mt(i_site - 1) == ID
                    n_hits_double = n_hits_double + 1;
                else
                    n_hits_single = n_hits_single + 1;
                end
            elseif (mt(i_site + 1) == ID || mt(i_site - 1) == ID)
                n_hits_double = n_hits_double + 1;
            else
                n_hits_single = n_hits_single + 1;
            end
        end
    end
end

n_hits_double = n_hits_double / 2;

single_fract = n_hits_single / (n_hits_double + n_hits_single);
double_fract = n_hits_double / (n_hits_double + n_hits_single);

data = [single_fract double_fract];
labels = {'Singly-bound', 'Doubly-bound'};
pie(data, labels);
ylabel('Fractional time in state')
title('Singly-bound vs doubly-bound kinesin')
%{
msg1 = sprintf('Fraction single-bound: %f', single_fract);
msg2 = sprintf('Fraction double-bound: %f', double_fract);
disp(msg1)
disp(msg2)
%}