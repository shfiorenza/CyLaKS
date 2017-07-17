clear all
n_timesteps = 100000;
length_of_microtubule = 160;

raw_data = zeros([length_of_microtubule 2*n_timesteps]);
temp_one = zeros([length_of_microtubule 1]);
mt_one = zeros([length_of_microtubule 1]);
final_data = zeros([length_of_microtubule 8]);

data_file = fopen('/home/shane/Desktop/code_cpp/test sims/60mil_160_LHLH.file');
raw_data = fread(data_file, [length_of_microtubule, 2*n_timesteps], '*int');
fclose(data_file);
raw_data((raw_data ~= 2) | (raw_data == 0)) = 0;
raw_data((raw_data == 2) | (raw_data ~= 0)) = 1;
for i=1:2:((2*n_timesteps)-1)
    temp_one(:, 1) = raw_data(:, i);
    mt_one(:, 1) = mt_one(:, 1) + double(temp_one(:, 1)./n_timesteps);
end
final_data(:, 1) = mt_one(:, 1);

clear raw_data;
clear mt_one; 
mt_one = zeros([length_of_microtubule 1]);
data_file = fopen('/home/shane/Desktop/code_cpp/test sims/60mil_160_Hn.file');
raw_data = fread(data_file, [length_of_microtubule, 2*n_timesteps], '*int');
fclose(data_file);
raw_data((raw_data ~= 2) | (raw_data == 0)) = 0;
raw_data((raw_data == 2) | (raw_data ~= 0)) = 1;
for i=1:2:((2*n_timesteps)-1)
    temp_one(:, 1) = raw_data(:, i);
    mt_one(:, 1) = mt_one(:, 1) + double(temp_one(:, 1)./n_timesteps);
end
final_data(:, 2) = mt_one(:, 1);

clear raw_data;
clear mt_one; 
mt_one = zeros([length_of_microtubule 1]);
data_file = fopen('/home/shane/Desktop/code_cpp/test sims/60mil_160_Hx.file');
raw_data = fread(data_file, [length_of_microtubule, 2*n_timesteps], '*int');
fclose(data_file);
raw_data((raw_data ~= 2) | (raw_data == 0)) = 0;
raw_data((raw_data == 2) | (raw_data ~= 0)) = 1;
for i=1:2:((2*n_timesteps)-1)
    temp_one(:, 1) = raw_data(:, i);
    mt_one(:, 1) = mt_one(:, 1) + double(temp_one(:, 1)./n_timesteps);
end
final_data(:, 3) = mt_one(:, 1);

clear raw_data;
clear mt_one; 
mt_one = zeros([length_of_microtubule 1]);
data_file = fopen('/home/shane/Desktop/code_cpp/test sims/60mil_160_LHx.file');
raw_data = fread(data_file, [length_of_microtubule, 2*n_timesteps], '*int');
fclose(data_file);
raw_data((raw_data ~= 2) | (raw_data == 0)) = 0;
raw_data((raw_data == 2) | (raw_data ~= 0)) = 1;
for i=1:2:((2*n_timesteps)-1)
    temp_one(:, 1) = raw_data(:, i);
    mt_one(:, 1) = mt_one(:, 1) + double(temp_one(:, 1)./n_timesteps);
end
final_data(:, 4) = mt_one(:, 1);

clear raw_data;
clear mt_one; 
mt_one = zeros([length_of_microtubule 1]);
data_file = fopen('/home/shane/Desktop/code_cpp/test sims/60mil_160_M.file');
raw_data = fread(data_file, [length_of_microtubule, 2*n_timesteps], '*int');
fclose(data_file);
raw_data((raw_data ~= 2) | (raw_data == 0)) = 0;
raw_data((raw_data == 2) | (raw_data ~= 0)) = 1;
for i=1:2:((2*n_timesteps)-1)
    temp_one(:, 1) = raw_data(:, i);
    mt_one(:, 1) = mt_one(:, 1) + double(temp_one(:, 1)./n_timesteps);
end
final_data(:, 5) = mt_one(:, 1);

clear raw_data;
clear mt_one; 
mt_one = zeros([length_of_microtubule 1]);
data_file = fopen('/home/shane/Desktop/code_cpp/test sims/60mil_160_Ln.file');
raw_data = fread(data_file, [length_of_microtubule, 2*n_timesteps], '*int');
fclose(data_file);
raw_data((raw_data ~= 2) | (raw_data == 0)) = 0;
raw_data((raw_data == 2) | (raw_data ~= 0)) = 1;
for i=1:2:((2*n_timesteps)-1)
    temp_one(:, 1) = raw_data(:, i);
    mt_one(:, 1) = mt_one(:, 1) + double(temp_one(:, 1)./n_timesteps);
end
final_data(:, 6) = mt_one(:, 1);

clear raw_data;
clear mt_one; 
mt_one = zeros([length_of_microtubule 1]);
data_file = fopen('/home/shane/Desktop/code_cpp/test sims/60mil_160_LHn.file');
raw_data = fread(data_file, [length_of_microtubule, 2*n_timesteps], '*int');
fclose(data_file);
raw_data((raw_data ~= 2) | (raw_data == 0)) = 0;
raw_data((raw_data == 2) | (raw_data ~= 0)) = 1;
for i=1:2:((2*n_timesteps)-1)
    temp_one(:, 1) = raw_data(:, i);
    mt_one(:, 1) = mt_one(:, 1) + double(temp_one(:, 1)./n_timesteps);
end
final_data(:, 7) = mt_one(:, 1);

clear raw_data;
clear mt_one; 
mt_one = zeros([length_of_microtubule 1]);
data_file = fopen('/home/shane/Desktop/code_cpp/test sims/60mil_160_Lx.file');
raw_data = fread(data_file, [length_of_microtubule, 2*n_timesteps], '*int');
fclose(data_file);
raw_data((raw_data ~= 2) | (raw_data == 0)) = 0;
raw_data((raw_data == 2) | (raw_data ~= 0)) = 1;
for i=1:2:((2*n_timesteps)-1)
    temp_one(:, 1) = raw_data(:, i);
    mt_one(:, 1) = mt_one(:, 1) + double(temp_one(:, 1)./n_timesteps);
end
final_data(:, 8) = mt_one(:, 1);

fig1 = figure(1);
set(fig1,'Position', [50, 50, 1.5*560, 1.5*420])
plot(linspace(0, length_of_microtubule, length_of_microtubule), final_data);
%style stuff
ylim([0 1])
grid on
grid minor
xlabel('Microtubule site')
ylabel('Fraction of the time occupied')