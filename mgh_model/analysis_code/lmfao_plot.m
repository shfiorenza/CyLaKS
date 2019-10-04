concentrations = [20,80,120,220,420];
runlengths = [2.0, 2.4, 2.5, 2.4 , 2.4];
lifetimes = [9.7, 11.4, 11.6, 11.3, 11.5];
velocities = [200, 210, 220, 210, 210];

exp_runlengths = [1.3, 2.8, 2.2, 3.0, 2.9];
exp_lifetimes = [4.0, 13.7, 14., 19.2, 26.0];
exp_velocities = [385, 289, 289, 293, 209];
    
subplot(1, 3, 1)
plot(concentrations, runlengths, '--o','LineWidth', 2);
hold on
plot(concentrations, exp_runlengths, ':d','LineWidth', 2);
ylabel('Run length (microns)');
ylim([0 4]);
%xlabel('KIF4A concentration (pM)');
subplot(1, 3, 2)
plot(concentrations, lifetimes, '--o', 'LineWidth', 2);
hold on
plot(concentrations, exp_lifetimes, ':d','LineWidth', 2);
ylabel('Lifetime (seconds)');
ylim([0 30]);
xlabel('KIF4A concentration (pM)');
subplot(1, 3, 3)
plot(concentrations, velocities, '--o', 'LineWidth', 2);
hold on
plot(concentrations, exp_velocities, ':d','LineWidth', 2);
ylabel('Velocity (nm/s)');
ylim([0 500]);
%xlabel('KIF4A concentration (pM)');
