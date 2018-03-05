clear all
n_datapoints = 100000;
n_mt_lengths = 4;
species_ID = 2;

fileDirectory = '/home/shane/Projects/overlap_analysis/mgh_model/%s';
nameStruct = '%#.2f_%#.1f_%i_occupancy.file';

% Run through different failstep rates
for i=1:1:1
   failstep_rate = i * 2.5;
   % For each failstep rate, run through different k_off_pseudo values
   endtag_data = [(n_mt_lengths) 10];
   for j=1:1:10
       k_off_pseudo = j / 2;
       % Scan through desired MT lengths (2um, 4um, 6um, 10um)
       for k = 1:1:5
           if(k ~= 4)
               mt_length = k * 250;
              
               temp_mt = zeros([mt_length 1]);
               final_mt = zeros([mt_length 1]);
             
               fileName = sprintf(nameStruct, k_off_pseudo, ...
                     failstep_rate, mt_length)
             
               data_file = fopen(sprintf(fileDirectory, fileName));
               raw_data = fread(data_file, [mt_length, n_datapoints], '*int');
               fclose(data_file);
              
               raw_data((raw_data ~= species_ID) | (raw_data == 0)) = 0;
               raw_data((raw_data == species_ID) | (raw_data ~= 0)) = 1;
            
               % Read in and average occupancy data over all datapoints
               for i=1:1:n_datapoints
                   temp_mt(:, 1) = raw_data(:, i);
                   final_mt(:, 1) = final_mt(:, 1) + double(temp_mt(:, 1)./n_datapoints);
               end
           
               % Find endtag length based on final MT occupancy data
               for i=1:1:mt_length
                   site_occupancy = final_mt(i, 1);
                   if(site_occupancy > 0.5)
                       endtag_site =  i;
                       break;
                   end
               end
               endtag_length = (mt_length - endtag_site) * 0.008;
               if(k == 5)
                   endtag_data(5, j) = endtag_length;
                   avg = (endtag_data(5, j) + endtag_data(3, j))/ 2;
                   endtag_data(4, j) = avg;
               else
                   endtag_data(k, j) = endtag_length;
               end
           end
       end    
   end
   plot(linspace(2, 10, 5), endtag_data, 'LineWidth', 2);
   legend('0.5', '1.0', '1.5', '2.0', '2.5', '3.0', '3.5', '4.0', '4.5', '5.0');
   legend('Location', 'northeastoutside');
   title(sprintf('failstep rate = %g', failstep_rate));
   xlabel("Length of microtubule (microns)");
   ylabel("Length of endtag (microns)");

end