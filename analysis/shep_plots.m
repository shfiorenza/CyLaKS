
vel = [506.1, 241.0, 140.0, 55.0, 54.5];
n_proto = [1, 2, 3, 5, 8];

%set(0, 'DefaultAxesFontName', 'Arial');
%set(0, 'DefaultTextFontName', 'Arial');

fig = figure('Position', [50 50 720 600]);
set(gcf, 'DefaultAxesFontName', 'Arial');
set(gcf, 'DefaultTextFontName', 'Arial');
hold on
plot(n_proto, vel, '.', 'MarkerSize', 50)

ylabel("Shepherding velocity (nm/s)");
yticks([0 250 500])
xlabel("Protofilament number")
xlim([0 10]);
xticks([1 2 3 5 8]);

set(gca,'box','off')
set(gca, 'FontSize', 28);
%fontname("arial");
%set(gca, 'FontName', 'arial');
set(gca,'TickDir','out');
set(gca,'LineWidth',1,'TickLength',[0.025 0.025]);