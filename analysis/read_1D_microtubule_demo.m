clear all;
aa = 0;
conc = [500 2000 4000];
velocity = 2.*1000./8;
one_D_array1 = zeros([length(conc) 4002]);
one_D_array2 = zeros([length(conc) 4002]);
time = 100000;
size = 600;
t_freq = 30;

for jj = conc;
    aa = aa + 1;
    formatSpec = '%01d';
    str = num2str(5,formatSpec);
    fid1= fopen(['/Users/hskuan/Documents/Results/two_seg/600_f_0_v_0_50_bc_', num2str(jj),'_', num2str(3300), '_five_seg_3300/detail.file_2']);
        
    data1 = zeros([1 4002]);
    data2 = zeros([1 4002]);
    
    for ii=1:1:time
        data1(1:1:4002) = fread(fid1, 4002, '*int');

        data2 = data1;
        
        data1((data1~=2) &(data1~=4)) = 0;
        data1((data1==2) | (data1==4)) = 1;
        data2((data2~=3) &(data2~=4)) = 0;
        data2((data2==3) | (data2==4)) = 1;
        if mod(ii, 1)==0
        one_D_array1(aa, :) = one_D_array1(aa, :) + double(data1(1:1:4002))./time;%0.5*(double(data1(1:1:4002)) + double(data2(1:1:4002)))./20000000;%
        
        one_D_array2(aa, :) = one_D_array2(aa, :) + double(data2(1:1:4002))./time;%0.5*(double(data1(1:1:4002)) - double(data2(1:1:4002)))./20000000;%
        end

    end
    
    fclose(fid1);
end

fig1 = figure(1);
set(fig1,'Position', [50, 50, 1.5*560, 1.5*420])
for ii=1:1:aa  %[1 2 4 7 10 13 16];%
plot(linspace(-size,size,2*size)./(2*size), ([one_D_array1(ii, (2000+size+1):-1:2002) one_D_array1(ii, 1:1:(size))]), '-.'); %right
    hold on;
plot(linspace(-size,size,2*size)./(2*size), ([one_D_array2(ii, (2000+size+1):-1:2002) one_D_array2(ii, 1:1:(size))]), '-.'); %left
    hold on;%./(2*size)
end
ylim([0 1])