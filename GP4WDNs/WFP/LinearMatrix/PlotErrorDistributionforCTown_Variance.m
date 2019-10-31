function PlotErrorDistributionforCTown_Variance(error)
fontsize = 80;
figure2 = figure;
% Create axes
%axes2 = axes('Parent',figure2);
h1=histogram(error,'Normalization','probability');

xlabel('$\mathrm{RE(\%)}$','FontSize',fontsize+3,'interpreter','latex');
ylabel('$\mathrm{Frequency}$','FontSize',fontsize+2,'interpreter','latex');
set(gca, 'TickLabelInterpreter', 'latex','fontsize',fontsize+10);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 24 8])
print(figure2,'Varianceerrordistribution','-depsc2','-r300');
end



