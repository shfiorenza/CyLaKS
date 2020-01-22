concentrations = [20, 50, 80, 120, 220, 420];

runlengths = [1.1, 1.0, 1.2, 1.2, 1.7, 3.9];
err_runlengths = [0.1, 0.07, 0.06, 0.05, 0.05, 0.09];
lifetimes = [2.3, 1.9, 2.9, 3.2, 5.1, 19.3];
err_lifetimes = [0.2, 0.1, 0.2, 0.1, 0.2, 0.4];
velocities = [490, 540, 420, 390, 330, 200];
err_velocities = [70, 50, 30, 20, 20, 10];

%{
runlengths = [1.1, 1.2, 1.6, 1.5, 1.8, 4.6];
lifetimes = [1.8, 2.4, 4.2, 3.9, 6.2, 27];
velocities = [610, 520, 380, 380, 300, 170];
%}

%{
runlengths = [1.0, 1.0, 1.1, 1.1, 1.5, 2.8];
lifetimes = [1.7, 1.8, 2.3, 3.4, 4.8, 13.2];
velocities = [560, 570, 470, 330, 320, 210];
%}

%{
runlengths(1, :) = [0.9, 1.4, 1.3,  1.3, 2.0, 4.3];
err_runlengths(1, :) = [0.1, 0.09, 0.07, 0.06, 0.07, 0.1];  
runlengths(2, :) = [0.9, 1.0, 1.2, 1.4, 1.9, 2.2];
err_runlengths(2, :) = [0.1, 0.07, 0.06, 0.07, 0.07, 0.09];

lifetimes(1, :) = [1.7, 2.3, 2.1, 2.3, 3.8, 8.3];  
err_lifeitmes(1, :) = [0.2, 0.2, 0.1, 0.1, 0.1, 0.2];
lifetimes(2, :) = [1.7, 3.7, 8.4, 16.0, 35.3, 84.5];
err_lifetimes(2, :) = [0.2, 0.3, 0.5, 0.7, 1, 3];

velocities(1, :) = [550, 610, 600, 600, 530, 520];
err_velocities(1, :) = [90, 60, 50, 40, 20, 20];
velocities(2, :) = [550, 280, 140, 90, 50, 30];
err_velocities(2, :) = [90, 30, 10, 10, 0, 0];
%}

exp_runlengths = [0.970, 1.310, 2.420, 1.660, 1.960, 2.86];
exp_err_runlengths = [0.180, 0.320, 0.350, 0.940, 0.310, 0.72];
exp_lifetimes = [1.8, 2.1, 7.1, 5.2, 8.3, 17.9];
exp_err_lifetimes = [0.6, 0.7, 1.7, 5.9, 2.6, 3.9];
exp_velocities = [600, 710, 360, 310, 310, 180];
exp_err_velocities = [75, 110, 50, 79, 40, 40];

fig1 = figure();
set(fig1, 'Position', [50, 50, 2*720, 2*240])
    
subplot(1, 4, 1)
errorbar(concentrations, exp_runlengths, exp_err_runlengths, ':d','LineWidth', 2);
hold on
errorbar(concentrations, runlengths, err_runlengths, '--o','LineWidth', 2);
ylabel('Run length (microns)');
xlim([0 440]);
ylim([0 5]);

subplot(1, 4, 2)
errorbar(concentrations, exp_lifetimes, exp_err_lifetimes, ':d','LineWidth', 2);
hold on
errorbar(concentrations, lifetimes, err_lifetimes, '--o', 'LineWidth', 2);
xlabel('KIF4A concentration (pM)');
ylabel('Life time (seconds)');
xlim([0 440]);
ylim([0 30]);

subplot(1, 4, 3)
errorbar(concentrations, exp_velocities, exp_err_velocities, ':d','LineWidth', 2);
hold on
errorbar(concentrations, velocities, err_velocities, '--o', 'LineWidth', 2);
ylabel('Velocity (nm/s)');
xlim([0 440]);
ylim([0 1000]);

subplot(1, 4, 4);
plot(0,0,  0,0,  0,0,  0,0)
axis off
legend({'Experimental data', 'Simulation data'}); %,  ...
    %'Sim data w/ cooperative unbind\_ii'}, 'location', 'best');
