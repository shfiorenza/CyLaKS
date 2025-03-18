
vel = [506.1, 241.0, 140.0, 55.0, 54.5];
n_proto = [1, 2, 3, 5, 8];

fig = figure('Position', [50 50 720 600]);
plot(n_proto, vel, '.', 'MarkerSize', 50)

set(gca,'box','off')
set(gca, 'FontSize', 28);
set(gca,'TickDir','out');
set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);

ylabel("Shepherding velocity (nm/s)");
yticks([0 250 500])
xlabel("Protofilament number")
xlim([0 10]);
xticks([1 2 3 5 8])