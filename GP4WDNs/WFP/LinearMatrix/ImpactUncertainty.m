close all
load('ImpactUncertainty.mat')
width = 0.3;
Figure1 = figure;
bartest = bar3(Bar3UncertaintyImpactCopyCopy);
hold on
xlim([0 28])
zlim([0 25])
%[1 3 5 8 10 12 14 16 18]
set(gca,'FontSize',25,'TickLabelInterpreter','latex','XTick', [1 4 7 11 14 17 21 24 27],...
    'XTickLabel',{'$\epsilon_{0\%}$','$\epsilon_{2.5\%}$','$\epsilon_{5\%}$','$d_{0\%}$','$d_{15\%}$','$d_{30\%}$','$c_{0\%}$','$c_{15\%}$','$c_{30\%}$'})
set(gca,'FontSize',25,'TickLabelInterpreter','latex','YTick', [4  8 13 17],...
    'YTickLabel',{'J2 $\sim$ J7','T8','P23 $\sim$ P78', 'PU12'})
xlabel('$\mathrm{Uncertainty}$','FontSize',25,'interpreter','latex');
%ylabel('$\mathrm{ID}$','FontSize',25,'interpreter','latex');
zlabel('$\mathrm{Standard\ deviation}:\sigma$','FontSize',25,'interpreter','latex');

pos = get(gca, 'Position');
pos(1) = 0.13;
pos(2) = 0.055;
pos(4) = 1.1;
set(gca, 'Position', pos)
view(gca,[-19.2 39.8]);