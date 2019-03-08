clear all;
n_sites = 1000;
initial_shift = [800,600,400,200,0];
%plus_end_dist = (-initial_shift + n_sites)*0.008; 
initial_overlap = (1000 - initial_shift)*0.008;

max_vel = [5,50,100,160,180];
%final_overlap_length = [1,1.1,1.2,1.2,1.4,1.4,1.4,1.4,1.3,1.4,0];


fig1 = figure();
set(fig1, 'Position', [50, 50, 640, 420]);
plot(initial_overlap, max_vel, '-o', 'LineWidth', 2);
title('Sliding for 8 micron-long microtubules', 'FontSize', 14);
ylabel('Maximum sliding velocity (nm/s)', 'FontName', 'Garuda','FontSize', 12);
xlabel('Initial overlap length (microns)', 'FontName', 'Garuda','FontSize', 12)
xticks([0,2,4,6,8]);
axes = gca;
axes.XAxis.FontSize = 14;
axes.YAxis.FontSize = 14;

%{
yyaxis right
plot(plus_end_dist, final_overlap_length, 'LineWidth', 2); 
ylabel('Final overlap length (microns)');
ylim([0 10]);
%}