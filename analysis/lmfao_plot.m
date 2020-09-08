concentrations = [20,80,120,220,420];
runlengths = [1.5, 1.5, 1.6, 1.5, 1.4];
lifetimes = [4.5, 4.8, 5.3, 6.6, 13.0];
velocities = [340, 320, 290, 230, 100];

subplot(1, 3, 1)
plot(concentrations, runlengths, '--o','LineWidth', 2);
ylabel('Run length (microns)');
ylim([0 4]);
%xlabel('KIF4A concentration (pM)');
subplot(1, 3, 2)
plot(concentrations, lifetimes, '--o', 'LineWidth', 2);
ylabel('Lifetime (seconds)');
ylim([0 30]);
xlabel('KIF4A concentration (pM)');
subplot(1, 3, 3)
plot(concentrations, velocities, '--o', 'LineWidth', 2);
ylabel('Velocity (nm/s)');
ylim([0 500]);
%xlabel('KIF4A concentration (pM)');
