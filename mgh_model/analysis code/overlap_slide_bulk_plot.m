clear all;
n_sites = 1000;
initial_shift = [-1000,-800,-600,-400,-200,0,200,400,600,800,1000];
plus_end_dist = (-initial_shift + n_sites)*0.008; 

max_vel = [700,300,275,175,225,160,115,90,60,35,0];
final_overlap_length = [1,1.1,1.2,1.2,1.4,1.4,1.4,1.4,1.3,1.4,0];


fig1 = figure();
set(fig1, 'Position', [50, 50, 960, 600]);

title('8-micron microtubules; varying initial overlap length');
yyaxis left
plot(plus_end_dist, max_vel, 'LineWidth', 2);
ylabel('Maximum sliding velocity (nm/s)');
xlabel('Initial distance between plus-ends of MTs (microns)')

yyaxis right
plot(plus_end_dist, final_overlap_length, 'LineWidth', 2); 
ylabel('Final overlap length (microns)');
ylim([0 10]);